#!/usr/bin/env python3
"""
Query Copernicus Data Space Ecosystem STAC for Sentinel-2 L2A observations at a
point and sample the Scene Classification Layer (SCL) at that point.

The script does not use Google Earth Engine.

Notes on "exact pixel acquisition time":
CDSE STAC exposes the observation/item acquisition timestamp. True detector-line
pixel timing is not present in the STAC item or SCL raster asset, so this script
reports the best available acquisition time from item metadata and labels the
source explicitly. If you need sub-scene detector timing, you must add a SAFE
metadata timing model for the specific band/detector; the SCL cloud class alone
does not contain that timing information.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import html
import http.server
import json
import math
import os
import re
import secrets
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import webbrowser
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Any, Iterable


STAC_SEARCH_URL = "https://stac.dataspace.copernicus.eu/v1/search"
CDSE_TOKEN_URL = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
DEFAULT_COLLECTION = "sentinel-2-l2a"

SCL_CLASSES = {
    0: "No data",
    1: "Saturated or defective",
    2: "Topographic cast shadows",
    3: "Cloud shadows",
    4: "Vegetation",
    5: "Not vegetated",
    6: "Water",
    7: "Unclassified",
    8: "Cloud medium probability",
    9: "Cloud high probability",
    10: "Thin cirrus",
    11: "Snow or ice",
}

CLOUD_SCL_VALUES = {8, 9, 10}
PROBLEMATIC_SCL_VALUES = {3, 8, 9, 10, 11}
TRAILING_PARITY_BY_SENSOR = {
    "S2A": "odd",
    "S2B": "odd",
    "S2C": "even",
}


@dataclass
class PixelSample:
    value: int | None
    row: int | None
    col: int | None
    asset_key: str | None
    asset_href: str | None
    error: str | None = None


@dataclass
class DetectorSample:
    detector_id: str | None
    source: str | None
    row: int | None = None
    col: int | None = None
    href: str | None = None
    error: str | None = None


@dataclass
class PixelTiming:
    time_utc: str | None
    source: str
    granule_sensing_time_utc: str | None = None
    granule_sensing_time_position: str | None = None
    detector_id: str | None = None
    detector_source: str | None = None
    detector_parity: str | None = None
    detector_trailing_parity: str | None = None
    detector_parallax_seconds: float | None = None
    band_id: str | None = None
    method: str | None = None
    candidate_count: int = 0
    candidate_min_utc: str | None = None
    candidate_max_utc: str | None = None
    candidate_times_json: str | None = None
    error: str | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "List Sentinel-2 L2A CDSE observations intersecting a latitude/"
            "longitude and sample pixel-wise cloud class from SCL."
        )
    )
    parser.add_argument("--lat", type=float, required=True, help="Latitude in EPSG:4326.")
    parser.add_argument("--lon", type=float, required=True, help="Longitude in EPSG:4326.")
    parser.add_argument(
        "--start",
        help="Optional UTC start datetime/date, e.g. 2024-01-01 or 2024-01-01T00:00:00Z.",
    )
    parser.add_argument(
        "--end",
        help="Optional UTC end datetime/date, e.g. 2024-12-31 or 2024-12-31T23:59:59Z.",
    )
    parser.add_argument(
        "--collection",
        default=DEFAULT_COLLECTION,
        help=f"STAC collection to query. Default: {DEFAULT_COLLECTION}.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=100,
        help="STAC page size. CDSE may cap this. Default: 100.",
    )
    parser.add_argument(
        "--max-items",
        type=int,
        default=0,
        help="Maximum items to process; 0 means no explicit limit. Default: 0.",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.0,
        help="Seconds to sleep between STAC pages. Default: 0.",
    )
    parser.add_argument(
        "--no-sample",
        action="store_true",
        help="Only list observations; do not open SCL rasters.",
    )
    parser.add_argument(
        "--band-id",
        default="4",
        help=(
            "Sentinel-2 band id used for detector datation. Default: 4 (B04). "
            "Use the numeric metadata band id, not B04 text."
        ),
    )
    parser.add_argument(
        "--detector",
        help=(
            "Optional detector id to force for acquisition-time extrapolation. "
            "If omitted, the script tries to sample the granule's CDSE MSK_DETFOO mask."
        ),
    )
    parser.add_argument("--detector-mask-s2a", help="Optional local S2A detector footprint mask JP2/TIF.")
    parser.add_argument("--detector-mask-s2b", help="Optional local S2B detector footprint mask JP2/TIF.")
    parser.add_argument("--detector-mask-s2c", help="Optional local S2C detector footprint mask JP2/TIF.")
    parser.add_argument(
        "--timing-origin",
        choices=("granule", "metadata", "scene-start", "scene-center"),
        default="granule",
        help=(
            "Pixel timing origin. 'granule' uses the studied tile's MTD_TL.xml SENSING_TIME. "
            "'scene-start' and 'scene-center' use STAC datetime. 'metadata' uses GPS_TIME/"
            "REFERENCE_LINE from SAFE datastrip metadata. Default: granule."
        ),
    )
    parser.add_argument(
        "--granule-sensing-time-position",
        choices=("center", "start"),
        default="start",
        help=(
            "How to interpret MTD_TL.xml SENSING_TIME in row-line timing. "
            "The Sentinel-2 PSD documents this as first-line/sensing-start time for Granule/Tile PDI. "
            "Default: start."
        ),
    )
    parser.add_argument(
        "--time-per-line",
        type=float,
        default=0.001515,
        help="Seconds per 10 m full-resolution Sentinel-2 line. Default: 0.001515.",
    )
    parser.add_argument(
        "--detector-parallax",
        type=float,
        default=2.6,
        help="Optional detector focal-plane parallax delay in seconds. Default: 2.6.",
    )
    parser.add_argument(
        "--trailing-detectors",
        choices=("odd", "even", "none"),
        default="none",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--metadata-root",
        help=(
            "Optional local directory containing downloaded .SAFE products. "
            "The script will use local MTD_DS.xml/MTD_TL.xml instead of remote assets."
        ),
    )
    parser.add_argument(
        "--cdse-token",
        default=os.environ.get("CDSE_ACCESS_TOKEN"),
        help="Optional CDSE bearer token. Defaults to CDSE_ACCESS_TOKEN environment variable.",
    )
    parser.add_argument(
        "--cdse-username",
        default=os.environ.get("CDSE_USERNAME"),
        help="Optional CDSE username. Defaults to CDSE_USERNAME environment variable.",
    )
    parser.add_argument(
        "--cdse-password",
        default=os.environ.get("CDSE_PASSWORD"),
        help="Optional CDSE password. Defaults to CDSE_PASSWORD environment variable.",
    )
    parser.add_argument(
        "--no-auth-prompt",
        action="store_true",
        help="Do not open a browser page to ask for a CDSE token when credentials are missing.",
    )
    parser.add_argument(
        "--auth-prompt-timeout",
        type=int,
        default=300,
        help="Seconds to wait for the browser token prompt. Default: 300.",
    )
    parser.add_argument(
        "--cache-dir",
        default=".cdse_cache",
        help=(
            "Directory for downloaded CDSE assets used when the server does not support "
            "GDAL HTTP range reads. Default: .cdse_cache."
        ),
    )
    parser.add_argument(
        "--out",
        default="sentinel2_pixel_cloud_results.csv",
        help="Output CSV path. Default: sentinel2_pixel_cloud_results.csv.",
    )
    parser.add_argument(
        "--jsonl",
        help="Optional JSON Lines output path with full per-item records.",
    )
    return parser.parse_args()


def normalize_datetime(value: str | None, end: bool = False) -> str | None:
    if not value:
        return None
    if "T" not in value:
        return f"{value}T{'23:59:59' if end else '00:00:00'}Z"
    if value.endswith("Z") or "+" in value[10:] or "-" in value[10:]:
        return value
    return value + "Z"


def stac_datetime(start: str | None, end: str | None) -> str | None:
    start_norm = normalize_datetime(start, end=False)
    end_norm = normalize_datetime(end, end=True)
    if start_norm and end_norm:
        return f"{start_norm}/{end_norm}"
    if start_norm:
        return f"{start_norm}/.."
    if end_norm:
        return f"../{end_norm}"
    return None


def http_json(url: str, payload: dict[str, Any] | None = None) -> dict[str, Any]:
    data = None
    headers = {"Accept": "application/json"}
    if payload is not None:
        data = json.dumps(payload).encode("utf-8")
        headers["Content-Type"] = "application/json"
    request = urllib.request.Request(url, data=data, headers=headers, method="POST" if data else "GET")
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"HTTP {exc.code} from {url}: {body[:1000]}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Could not reach {url}: {exc}") from exc


def get_cdse_token(username: str | None, password: str | None) -> str | None:
    if not username or not password:
        return None
    payload = urllib.parse.urlencode(
        {
            "client_id": "cdse-public",
            "username": username,
            "password": password,
            "grant_type": "password",
        }
    ).encode("utf-8")
    request = urllib.request.Request(
        CDSE_TOKEN_URL,
        data=payload,
        headers={"Accept": "application/json", "Content-Type": "application/x-www-form-urlencoded"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            token_response = json.loads(response.read().decode("utf-8"))
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"Could not obtain CDSE access token: {exc}") from exc
    token = token_response.get("access_token")
    if not token:
        raise RuntimeError("CDSE token response did not contain access_token.")
    return token


def extract_access_token(value: str) -> str | None:
    text = value.strip()
    if not text:
        return None
    if text.lower().startswith("bearer "):
        text = text[7:].strip()
    if text.startswith("{"):
        try:
            parsed = json.loads(text)
        except json.JSONDecodeError:
            return None
        token = parsed.get("access_token")
        return str(token).strip() if token else None
    return text


def ask_cdse_token_via_web_page(timeout_seconds: int = 300) -> str | None:
    token_result: dict[str, str] = {}
    ready = threading.Event()
    nonce = secrets.token_urlsafe(18)

    class TokenPromptHandler(http.server.BaseHTTPRequestHandler):
        def log_message(self, format: str, *args: Any) -> None:  # noqa: A002
            return

        def send_text(self, status: int, body: str, content_type: str = "text/html; charset=utf-8") -> None:
            encoded = body.encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(encoded)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(encoded)

        def do_GET(self) -> None:  # noqa: N802
            if self.path != f"/{nonce}":
                self.send_text(404, "Not found", "text/plain; charset=utf-8")
                return
            escaped_token_url = html.escape(CDSE_TOKEN_URL)
            body = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>CDSE Access Token</title>
  <style>
    body {{ font-family: system-ui, sans-serif; max-width: 760px; margin: 48px auto; line-height: 1.45; }}
    textarea {{ box-sizing: border-box; width: 100%; min-height: 180px; font-family: ui-monospace, monospace; }}
    button {{ padding: 10px 16px; margin-top: 12px; }}
    code {{ background: #f2f2f2; padding: 2px 4px; }}
  </style>
</head>
<body>
  <h1>Paste CDSE Access Token</h1>
  <p>The script needs a Copernicus Data Space Ecosystem bearer token to read OData assets.</p>
  <p>Paste either the raw token, <code>Bearer &lt;token&gt;</code>, or a JSON response containing <code>access_token</code>.</p>
  <form method="post" action="/{nonce}">
    <textarea name="token" autocomplete="off" spellcheck="false" autofocus></textarea>
    <br>
    <button type="submit">Use Token</button>
  </form>
  <p>Token endpoint used by CDSE tools: <code>{escaped_token_url}</code></p>
  <p>You can close this browser tab after submitting.</p>
</body>
</html>"""
            self.send_text(200, body)

        def do_POST(self) -> None:  # noqa: N802
            if self.path != f"/{nonce}":
                self.send_text(404, "Not found", "text/plain; charset=utf-8")
                return
            length = int(self.headers.get("Content-Length", "0"))
            payload = self.rfile.read(length).decode("utf-8", errors="replace")
            fields = urllib.parse.parse_qs(payload)
            token = extract_access_token(fields.get("token", [""])[0])
            if not token:
                self.send_text(400, "<p>No access token was provided.</p>")
                return
            token_result["token"] = token
            ready.set()
            self.send_text(200, "<p>Token received. You can close this tab and return to the script.</p>")

    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), TokenPromptHandler)
    url = f"http://127.0.0.1:{server.server_port}/{nonce}"
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    print(f"No CDSE credentials found. Opening token prompt: {url}", file=sys.stderr)
    try:
        webbrowser.open(url)
    except Exception as exc:  # noqa: BLE001
        print(f"Could not open browser automatically: {exc}", file=sys.stderr)
    try:
        ready.wait(timeout_seconds)
    finally:
        server.shutdown()
        server.server_close()
    return token_result.get("token")


def needs_cdse_token(args: argparse.Namespace) -> bool:
    return not args.no_sample


def resolve_cdse_token(args: argparse.Namespace) -> str | None:
    if args.cdse_token:
        return extract_access_token(args.cdse_token)
    token = get_cdse_token(args.cdse_username, args.cdse_password)
    if token:
        return token
    if needs_cdse_token(args) and not args.no_auth_prompt:
        return ask_cdse_token_via_web_page(args.auth_prompt_timeout)
    return None


def http_bytes(url: str, token: str | None = None) -> bytes:
    headers = {"Accept": "application/xml,text/xml,*/*"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    request = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            return response.read()
    except urllib.error.HTTPError as exc:
        body = exc.read(500).decode("utf-8", errors="replace")
        raise RuntimeError(f"HTTP {exc.code} from {url}: {body}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Could not reach {url}: {exc}") from exc


def safe_filename(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def asset_extension(href: str, default: str = ".bin") -> str:
    path = urllib.parse.urlparse(href).path
    base = os.path.basename(path)
    for suffix in (".jp2", ".tif", ".tiff", ".xml", ".json", ".zip"):
        if base.lower().endswith(suffix):
            return suffix
    return default


def download_asset_to_cache(
    item: dict[str, Any],
    asset_key: str,
    href: str,
    token: str | None,
    cache_dir: str,
    default_extension: str = ".bin",
) -> str:
    os.makedirs(cache_dir, exist_ok=True)
    filename = safe_filename(f"{item.get('id', 'item')}_{asset_key}") + asset_extension(href, default_extension)
    path = os.path.join(cache_dir, filename)
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path

    headers = {"Accept": "*/*"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    request = urllib.request.Request(href, headers=headers)
    tmp_path = path + ".part"
    try:
        with urllib.request.urlopen(request, timeout=300) as response, open(tmp_path, "wb") as handle:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
        os.replace(tmp_path, path)
    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise
    return path


def iter_stac_items(args: argparse.Namespace) -> Iterable[dict[str, Any]]:
    payload: dict[str, Any] = {
        "collections": [args.collection],
        "intersects": {"type": "Point", "coordinates": [args.lon, args.lat]},
        "limit": args.limit,
        "sortby": [{"field": "datetime", "direction": "asc"}],
    }
    dt = stac_datetime(args.start, args.end)
    if dt:
        payload["datetime"] = dt

    url = STAC_SEARCH_URL
    yielded = 0

    while url:
        page = http_json(url, payload)
        payload = None

        for item in page.get("features", []):
            yield item
            yielded += 1
            if args.max_items and yielded >= args.max_items:
                return

        next_url = None
        next_payload = None
        for link in page.get("links", []):
            if link.get("rel") == "next":
                next_url = link.get("href")
                if str(link.get("method", "")).upper() == "POST":
                    next_payload = link.get("body")
                break
        url = next_url
        payload = next_payload
        if url and args.sleep:
            time.sleep(args.sleep)


def find_scl_asset(item: dict[str, Any]) -> tuple[str | None, str | None]:
    assets = item.get("assets", {})
    preferred_keys = ("SCL_20m", "SCL", "scl", "SCL_60m", "scene_classification", "classification")
    for key in preferred_keys:
        href = asset_href(item, key)
        if href:
            return key, href
    for key, asset in assets.items():
        title = str(asset.get("title", "")).lower()
        roles = {str(role).lower() for role in asset.get("roles", [])}
        href = asset.get("href")
        if href and ("scl" in key.lower() or "scene classification" in title or "classification" in roles):
            return key, href
    return None, None


def asset_href(item: dict[str, Any], key: str) -> str | None:
    asset = item.get("assets", {}).get(key, {})
    alternate = asset.get("alternate", {}).get("https", {})
    return alternate.get("href") or asset.get("href")


def asset_requires_auth(item: dict[str, Any], key: str) -> bool:
    asset = item.get("assets", {}).get(key, {})
    refs = set(asset.get("auth:refs", []))
    alternate = asset.get("alternate", {}).get("https", {})
    refs.update(alternate.get("auth:refs", []))
    return bool(refs)


def asset_local_path(item: dict[str, Any], key: str) -> str | None:
    return item.get("assets", {}).get(key, {}).get("file:local_path")


def product_node_href(item: dict[str, Any], relative_parts: list[str]) -> str | None:
    product_href = item.get("assets", {}).get("Product", {}).get("href")
    if not product_href:
        return None
    base = product_href.rsplit("/$value", 1)[0]
    for part in relative_parts:
        base += f"/Nodes({part})"
    return base + "/$value"


def granule_relative_parts(item: dict[str, Any]) -> list[str] | None:
    local_path = asset_local_path(item, "granule_metadata")
    if not local_path:
        return None
    parts = local_path.replace("\\", "/").split("/")
    try:
        granule_index = parts.index("GRANULE")
    except ValueError:
        return None
    return parts[granule_index:]


def detector_mask_band_name(band_id: str) -> str:
    normalized = norm_metadata_id(band_id) or str(band_id)
    try:
        return f"B{int(normalized):02d}"
    except ValueError:
        text = str(band_id).upper()
        return text if text.startswith("B") else f"B{text}"


def find_local_metadata(metadata_root: str, item: dict[str, Any], key: str) -> str | None:
    local_path = asset_local_path(item, key)
    if not local_path:
        return None
    product_name = local_path.split(".SAFE", 1)[0] + ".SAFE"
    candidates = [
        os.path.join(metadata_root, local_path),
        os.path.join(metadata_root, product_name, local_path.split(".SAFE", 1)[-1].lstrip("/\\")),
    ]
    for candidate in candidates:
        if os.path.exists(candidate):
            return candidate
    for root, _dirs, files in os.walk(metadata_root):
        target = "MTD_DS.xml" if key == "datastrip_metadata" else "MTD_TL.xml"
        if target in files and product_name in root:
            return os.path.join(root, target)
    return None


def read_metadata_xml(
    item: dict[str, Any],
    key: str,
    metadata_root: str | None,
    token: str | None,
) -> ET.Element | None:
    if metadata_root:
        local = find_local_metadata(metadata_root, item, key)
        if local:
            return ET.parse(local).getroot()

    href = asset_href(item, key)
    if not href or href.startswith("s3://"):
        return None
    if asset_requires_auth(item, key) and not token:
        return None
    return ET.fromstring(http_bytes(href, token))


def granule_sensing_time(item: dict[str, Any], metadata_root: str | None, token: str | None) -> str | None:
    try:
        tl_root = read_metadata_xml(item, "granule_metadata", metadata_root, token)
    except Exception:
        return None
    if tl_root is None:
        return None
    return first_descendant_text(tl_root, {"SENSING_TIME", "Sensing_Time", "sensing_time"})


def lname(element: ET.Element) -> str:
    return element.tag.rsplit("}", 1)[-1]


def text_of(element: ET.Element | None) -> str | None:
    if element is None or element.text is None:
        return None
    value = element.text.strip()
    return value or None


def first_descendant_text(root: ET.Element, names: set[str]) -> str | None:
    for element in root.iter():
        if lname(element) in names:
            value = text_of(element)
            if value:
                return value
    return None


def parse_float(value: str | None) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_datetime(value: str | None) -> dt.datetime | None:
    if value is None:
        return None
    cleaned = value.strip().replace("Z", "+00:00")
    try:
        parsed = dt.datetime.fromisoformat(cleaned)
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=dt.timezone.utc)
    return parsed.astimezone(dt.timezone.utc)


def isoformat_utc(value: dt.datetime | None) -> str | None:
    if value is None:
        return None
    return value.astimezone(dt.timezone.utc).isoformat(timespec="microseconds").replace("+00:00", "Z")


def detector_parity(detector_id: str | None) -> str | None:
    detector = norm_metadata_id(detector_id)
    if detector is None:
        return None
    try:
        is_odd = int(detector) % 2 != 0
    except ValueError:
        return None
    return "odd" if is_odd else "even"


def automatic_detector_parallax(
    item: dict[str, Any],
    detector_id: str | None,
    parallax_seconds: float,
) -> tuple[float, str | None, str | None]:
    parity = detector_parity(detector_id)
    sensor_id = sensor_id_from_item(item)
    trailing_parity = TRAILING_PARITY_BY_SENSOR.get(sensor_id or "")
    if parity is None or trailing_parity is None:
        return 0.0, parity, trailing_parity
    return (parallax_seconds if parity == trailing_parity else 0.0), parity, trailing_parity


def child_texts(element: ET.Element) -> dict[str, str]:
    values: dict[str, str] = {}
    for child in list(element):
        value = text_of(child)
        if value:
            values[lname(child)] = value
    values.update({k.rsplit("}", 1)[-1]: v for k, v in element.attrib.items()})
    return values


def norm_metadata_id(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return str(int(text))
    except ValueError:
        return text


def nearest_metadata_value(
    element: ET.Element,
    parent_map: dict[ET.Element, ET.Element],
    names: tuple[str, ...],
) -> str | None:
    current: ET.Element | None = element
    while current is not None:
        values = child_texts(current)
        for name in names:
            if name in values:
                return values[name]
        current = parent_map.get(current)
    return None


def extract_timing_entries(root: ET.Element, band_id: str) -> list[dict[str, Any]]:
    parent_map = {child: parent for parent in root.iter() for child in list(parent)}
    wanted_band = norm_metadata_id(band_id)
    entries: list[dict[str, Any]] = []
    for detector in root.iter():
        values = child_texts(detector)
        names = {lname(child) for child in list(detector)}
        if not {"REFERENCE_LINE", "GPS_TIME"}.issubset(names):
            continue
        current_band = nearest_metadata_value(
            detector,
            parent_map,
            ("bandId", "band_id", "BAND_ID", "band", "band_id_1", "bandId_1"),
        )
        current_detector = (
            values.get("detectorId")
            or values.get("detector_id")
            or values.get("DETECTOR_ID")
            or values.get("detector")
            or values.get("Detector_Id")
        )
        current_band = norm_metadata_id(current_band)
        current_detector = norm_metadata_id(current_detector)
        if current_band is not None and current_band != wanted_band:
            continue
        entries.append(
            {
                "band_id": str(current_band or wanted_band or band_id),
                "detector_id": current_detector,
                "reference_line": parse_float(values.get("REFERENCE_LINE")),
                "gps_time": parse_datetime(values.get("GPS_TIME")),
            }
        )
    return [entry for entry in entries if entry["reference_line"] is not None and entry["gps_time"] is not None]


def extract_granule_positions(root: ET.Element, item: dict[str, Any]) -> dict[str, float]:
    local = asset_local_path(item, "granule_metadata") or ""
    granule_hint = os.path.basename(os.path.dirname(local.replace("\\", "/")))
    positions: dict[str, float] = {}
    for element in root.iter():
        values = child_texts(element)
        position = parse_float(values.get("POSITION"))
        if position is None:
            continue
        detector_id = (
            values.get("detectorId")
            or values.get("DETECTOR_ID")
            or values.get("detector_id")
            or values.get("detector")
        )
        text_blob = " ".join(value for value in values.values())
        if granule_hint and granule_hint not in text_blob:
            child_blob = " ".join((text_of(child) or "") for child in element.iter())
            if granule_hint not in child_blob:
                continue
        if detector_id is not None:
            positions[str(detector_id)] = position
    return positions


def estimate_pixel_timing(
    item: dict[str, Any],
    sample: PixelSample,
    band_id: str,
    detector_id: str | None,
    metadata_root: str | None,
    token: str | None,
    timing_origin: str,
    time_per_line: float,
    parallax_seconds: float,
    granule_sensing_position: str,
    detector_source: str | None,
) -> PixelTiming:
    properties = item.get("properties", {})
    fallback_time = first_property(properties, ("datetime", "start_datetime", "created"))
    granule_time = granule_sensing_time(item, metadata_root, token)
    if sample.row is None:
        return PixelTiming(
            None,
            "not computed",
            granule_sensing_time_utc=granule_time,
            granule_sensing_time_position=granule_sensing_position if granule_time else None,
            error=(
                "No sampled pixel row for timing extrapolation. "
                f"Scene acquisition time is available as acquisition_time_utc={fallback_time}."
            ),
        )

    # SCL_20m row to 10 m full-resolution line, using pixel centre.
    tile_line_10m = (sample.row + 0.5) * 2.0

    if timing_origin in {"granule", "scene-start", "scene-center"}:
        base_text = granule_time if timing_origin == "granule" else fallback_time
        base = parse_datetime(base_text)
        if base is None:
            return PixelTiming(
                None,
                "not computed",
                granule_sensing_time_utc=granule_time,
                granule_sensing_time_position=granule_sensing_position if granule_time else None,
                error=f"Could not parse {'granule SENSING_TIME' if timing_origin == 'granule' else 'STAC scene acquisition time'}.",
            )
        line_offset = tile_line_10m
        if timing_origin == "scene-center" or (timing_origin == "granule" and granule_sensing_position == "center"):
            line_offset -= 10980.0 / 2.0
        parallax, parity, trailing_parity = automatic_detector_parallax(item, detector_id, parallax_seconds)
        acq_time = base + dt.timedelta(seconds=line_offset * time_per_line + parallax)
        return PixelTiming(
            isoformat_utc(acq_time),
            f"{'granule MTD_TL.xml SENSING_TIME' if timing_origin == 'granule' else 'STAC ' + timing_origin} plus row-line model",
            granule_sensing_time_utc=granule_time,
            granule_sensing_time_position=granule_sensing_position if granule_time else None,
            detector_id=norm_metadata_id(detector_id),
            detector_source=detector_source,
            detector_parity=parity,
            detector_trailing_parity=trailing_parity,
            detector_parallax_seconds=parallax,
            band_id=band_id,
            method=(
                f"{base_text} + ((SCL_row + 0.5) * 2"
                f"{' - 5490' if timing_origin == 'scene-center' or (timing_origin == 'granule' and granule_sensing_position == 'center') else ''}) * {time_per_line}"
                f" + auto_detector_parallax(sensor={sensor_id_from_item(item)}, detector={norm_metadata_id(detector_id)}, "
                f"parity={parity}, trailing={trailing_parity}, seconds={parallax})"
            ),
            candidate_count=1,
            candidate_times_json=json.dumps(
                [
                    {
                        "detector_id": norm_metadata_id(detector_id) or "",
                        "detector_parity": parity or "",
                        "trailing_parity": trailing_parity or "",
                        "detector_parallax_seconds": parallax,
                        "band_id": str(band_id),
                        "time_utc": isoformat_utc(acq_time) or "",
                    }
                ],
                separators=(",", ":"),
            ),
        )

    try:
        ds_root = read_metadata_xml(item, "datastrip_metadata", metadata_root, token)
    except Exception as exc:  # noqa: BLE001
        return PixelTiming(
            None,
            "not computed",
            granule_sensing_time_utc=granule_time,
            granule_sensing_time_position=granule_sensing_position if granule_time else None,
            error=(
                f"Could not read datastrip metadata: {exc}. "
                f"Scene acquisition time is available as acquisition_time_utc={fallback_time}."
            ),
        )

    if ds_root is None:
        return PixelTiming(
            None,
            "not computed",
            granule_sensing_time_utc=granule_time,
            granule_sensing_time_position=granule_sensing_position if granule_time else None,
            error=(
                "Datastrip metadata is not locally available and the remote asset requires CDSE authentication. "
                f"Scene acquisition time is available as acquisition_time_utc={fallback_time}."
            ),
        )

    line_period = parse_float(first_descendant_text(ds_root, {"THEORETICAL_LINE_PERIOD", "LINE_PERIOD"}))
    if line_period is None:
        return PixelTiming(
            None,
            "not computed",
            granule_sensing_time_utc=granule_time,
            granule_sensing_time_position=granule_sensing_position if granule_time else None,
            error="No LINE_PERIOD or THEORETICAL_LINE_PERIOD in datastrip metadata.",
        )
    if line_period > 1:
        line_period /= 1000.0

    entries = extract_timing_entries(ds_root, band_id)
    if detector_id:
        normalized_detector_id = norm_metadata_id(detector_id)
        entries = [entry for entry in entries if norm_metadata_id(entry["detector_id"]) == normalized_detector_id]
    if not entries:
        return PixelTiming(
            None,
            "not computed",
            granule_sensing_time_utc=granule_time,
            granule_sensing_time_position=granule_sensing_position if granule_time else None,
            error=f"No band {band_id} detector timing entries found.",
        )

    candidates: list[dict[str, str]] = []
    for entry in entries:
        det = entry["detector_id"]
        datastrip_line = tile_line_10m
        delta_lines = datastrip_line - entry["reference_line"]
        parallax, parity, trailing_parity = automatic_detector_parallax(item, det, parallax_seconds)
        acq_time = entry["gps_time"] + dt.timedelta(seconds=delta_lines * line_period + parallax)
        candidates.append(
            {
                "detector_id": det or "",
                "detector_parity": parity or "",
                "trailing_parity": trailing_parity or "",
                "detector_parallax_seconds": str(parallax),
                "band_id": entry["band_id"],
                "time_utc": isoformat_utc(acq_time) or "",
            }
        )

    parsed_times = [parse_datetime(candidate["time_utc"]) for candidate in candidates]
    parsed_times = [value for value in parsed_times if value is not None]
    if not parsed_times:
        return PixelTiming(
            None,
            "not computed",
            granule_sensing_time_utc=granule_time,
            granule_sensing_time_position=granule_sensing_position if granule_time else None,
            error="Detector timing candidates could not be parsed.",
        )

    if len(candidates) == 1:
        return PixelTiming(
            candidates[0]["time_utc"],
            "SAFE datastrip detector datation model",
            granule_sensing_time_utc=granule_time,
            granule_sensing_time_position=granule_sensing_position if granule_time else None,
            detector_id=candidates[0]["detector_id"],
            detector_source=detector_source,
            detector_parity=candidates[0]["detector_parity"],
            detector_trailing_parity=candidates[0]["trailing_parity"],
            detector_parallax_seconds=float(candidates[0]["detector_parallax_seconds"]),
            band_id=candidates[0]["band_id"],
            method="GPS_TIME + (SCL_row_center_10m - REFERENCE_LINE) * LINE_PERIOD + detector_parallax",
            candidate_count=1,
            candidate_times_json=json.dumps(candidates, separators=(",", ":")),
        )

    return PixelTiming(
        None,
        "SAFE datastrip detector datation model; detector not resolved",
        granule_sensing_time_utc=granule_time,
        granule_sensing_time_position=granule_sensing_position if granule_time else None,
        band_id=band_id,
        method="Computed one candidate per available detector; pass --detector for a single extrapolated timestamp.",
        candidate_count=len(candidates),
        candidate_min_utc=isoformat_utc(min(parsed_times)),
        candidate_max_utc=isoformat_utc(max(parsed_times)),
        candidate_times_json=json.dumps(candidates, separators=(",", ":")),
    )


def sample_raster_pixel(rasterio_module: Any, transform_func: Any, href: str, lon: float, lat: float) -> tuple[int, int, int]:
    with rasterio_module.open(href) as dataset:
        xs, ys = transform_func("EPSG:4326", dataset.crs, [lon], [lat])
        row, col = dataset.index(xs[0], ys[0])
        if row < 0 or col < 0 or row >= dataset.height or col >= dataset.width:
            raise ValueError(f"Point is outside SCL raster bounds at row={row}, col={col}.")
        value = dataset.read(1, window=((row, row + 1), (col, col + 1)), boundless=False)[0, 0]
        return int(value), row, col


def sensor_id_from_item(item: dict[str, Any]) -> str | None:
    platform = str(item.get("properties", {}).get("platform", "")).upper()
    item_id = str(item.get("id", "")).upper()
    text = platform + " " + item_id
    if "SENTINEL-2A" in text or "S2A" in text:
        return "S2A"
    if "SENTINEL-2B" in text or "S2B" in text:
        return "S2B"
    if "SENTINEL-2C" in text or "S2C" in text:
        return "S2C"
    return None


def detector_mask_hrefs(item: dict[str, Any], band_id: str) -> list[tuple[str, str]]:
    granule_parts = granule_relative_parts(item)
    if not granule_parts:
        return []
    granule_dir = granule_parts[:-1]
    preferred_band = detector_mask_band_name(band_id)
    bands = [preferred_band] + [band for band in ("B04", "B02", "B03", "B08") if band != preferred_band]
    candidates: list[tuple[str, str]] = []
    for band in bands:
        relative_parts = granule_dir + ["QI_DATA", f"MSK_DETFOO_{band}.jp2"]
        href = product_node_href(item, relative_parts)
        if href:
            candidates.append((f"MSK_DETFOO_{band}", href))
    return candidates


def sample_detector_from_href(
    item: dict[str, Any],
    label: str,
    href: str,
    args: argparse.Namespace,
    token: str | None,
) -> DetectorSample:
    try:
        import rasterio
        from rasterio.warp import transform
    except ImportError:
        return DetectorSample(None, None, error="Install rasterio to sample detector masks.")

    if not token:
        return DetectorSample(None, None, href=href, error="Detector mask download requires CDSE authentication.")

    try:
        local_href = download_asset_to_cache(item, label, href, token, args.cache_dir, ".jp2")
        value, row, col = sample_raster_pixel(rasterio, transform, local_href, args.lon, args.lat)
        detector = norm_metadata_id(value)
        if detector in {None, "0"}:
            return DetectorSample(None, label, row=row, col=col, href=local_href, error="Detector mask sampled nodata/0.")
        return DetectorSample(detector, f"granule {label}", row=row, col=col, href=local_href)
    except Exception as exc:  # noqa: BLE001
        return DetectorSample(None, label, href=href, error=str(exc))


def resolve_detector_id(item: dict[str, Any], args: argparse.Namespace, token: str | None) -> DetectorSample:
    if args.detector:
        return DetectorSample(norm_metadata_id(args.detector), "forced --detector")

    last_error = None
    for label, href in detector_mask_hrefs(item, args.band_id):
        sample = sample_detector_from_href(item, label, href, args, token)
        if sample.detector_id:
            return sample
        last_error = sample.error

    sensor_id = sensor_id_from_item(item)
    mask_path = {
        "S2A": args.detector_mask_s2a,
        "S2B": args.detector_mask_s2b,
        "S2C": args.detector_mask_s2c,
    }.get(sensor_id or "")
    if not mask_path:
        return DetectorSample(None, None, error=last_error or "No detector mask sampled.")
    try:
        import rasterio
        from rasterio.warp import transform

        value, row, col = sample_raster_pixel(rasterio, transform, mask_path, args.lon, args.lat)
        return DetectorSample(norm_metadata_id(value), f"local {sensor_id} detector mask", row=row, col=col, href=mask_path)
    except Exception as exc:  # noqa: BLE001
        return DetectorSample(None, None, href=mask_path, error=f"Could not sample detector mask {mask_path}: {exc}")


def sample_scl(
    item: dict[str, Any],
    lon: float,
    lat: float,
    token: str | None = None,
    cache_dir: str = ".cdse_cache",
) -> PixelSample:
    key, href = find_scl_asset(item)
    if not href:
        return PixelSample(None, None, None, key, href, "No SCL asset found in STAC item.")
    if key and asset_requires_auth(item, key) and not token:
        return PixelSample(
            None,
            None,
            None,
            key,
            href,
            "SCL asset requires CDSE authentication. Set CDSE_ACCESS_TOKEN or CDSE_USERNAME/CDSE_PASSWORD.",
        )

    try:
        import rasterio
        from rasterio.warp import transform
    except ImportError as exc:
        return PixelSample(
            None,
            None,
            None,
            key,
            href,
            "Install rasterio to sample SCL rasters: python -m pip install rasterio",
        )

    try:
        env_options = {"GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR"}
        if token:
            env_options["GDAL_HTTP_HEADERS"] = f"Authorization: Bearer {token}"
        with rasterio.Env(**env_options):
            value, row, col = sample_raster_pixel(rasterio, transform, href, lon, lat)
            return PixelSample(value, row, col, key, href)
    except Exception as exc:  # noqa: BLE001 - preserve operational error details for users.
        message = str(exc)
        if "Range downloading not supported" not in message:
            return PixelSample(None, None, None, key, href, message)
        try:
            local_href = download_asset_to_cache(item, key or "SCL", href, token, cache_dir, ".jp2")
            value, row, col = sample_raster_pixel(rasterio, transform, local_href, lon, lat)
            return PixelSample(value, row, col, key, local_href)
        except Exception as retry_exc:  # noqa: BLE001
            return PixelSample(None, None, None, key, href, f"{message}; local-cache retry failed: {retry_exc}")


def first_property(properties: dict[str, Any], names: tuple[str, ...]) -> Any:
    for name in names:
        if name in properties:
            return properties[name]
    return None


def item_record(
    item: dict[str, Any],
    sample: PixelSample,
    timing: PixelTiming,
    lon: float,
    lat: float,
) -> dict[str, Any]:
    properties = item.get("properties", {})
    scl_value = sample.value
    scl_label = SCL_CLASSES.get(scl_value, "Unknown") if scl_value is not None else None
    item_cloud = first_property(properties, ("eo:cloud_cover", "cloudCover", "cloud_cover"))
    acquisition_time = first_property(properties, ("datetime", "start_datetime", "created"))
    acquisition_dt = parse_datetime(acquisition_time)
    granule_dt = parse_datetime(timing.granule_sensing_time_utc)
    granule_minus_stac_seconds = (
        (granule_dt - acquisition_dt).total_seconds() if acquisition_dt is not None and granule_dt is not None else None
    )

    return {
        "item_id": item.get("id"),
        "collection": item.get("collection"),
        "latitude": lat,
        "longitude": lon,
        "acquisition_time_utc": acquisition_time,
        "acquisition_time_source": "STAC item property; not detector-line pixel timing",
        "granule_sensing_time_utc": timing.granule_sensing_time_utc,
        "granule_sensing_time_position": timing.granule_sensing_time_position,
        "granule_minus_stac_seconds": granule_minus_stac_seconds,
        "pixel_acquisition_time_utc": timing.time_utc,
        "pixel_acquisition_time_source": timing.source,
        "pixel_acquisition_time_method": timing.method,
        "pixel_acquisition_time_band_id": timing.band_id,
        "pixel_acquisition_time_detector_id": timing.detector_id,
        "pixel_acquisition_time_detector_source": timing.detector_source,
        "pixel_acquisition_time_detector_parity": timing.detector_parity,
        "pixel_acquisition_time_detector_trailing_parity": timing.detector_trailing_parity,
        "pixel_acquisition_time_detector_parallax_seconds": timing.detector_parallax_seconds,
        "pixel_acquisition_time_candidate_count": timing.candidate_count,
        "pixel_acquisition_time_min_utc": timing.candidate_min_utc,
        "pixel_acquisition_time_max_utc": timing.candidate_max_utc,
        "pixel_acquisition_time_candidates_json": timing.candidate_times_json,
        "pixel_acquisition_time_error": timing.error,
        "platform": first_property(properties, ("platform", "constellation")),
        "product_type": first_property(properties, ("product:type", "s2:product_type")),
        "mgrs_tile": first_property(properties, ("grid:code", "s2:mgrs_tile")),
        "item_cloud_cover_percent": item_cloud,
        "scl_value": scl_value,
        "scl_label": scl_label,
        "pixel_is_cloud": scl_value in CLOUD_SCL_VALUES if scl_value is not None else None,
        "pixel_is_cloud_shadow": scl_value == 3 if scl_value is not None else None,
        "pixel_is_snow_or_ice": scl_value == 11 if scl_value is not None else None,
        "pixel_has_cloud_related_issue": scl_value in PROBLEMATIC_SCL_VALUES if scl_value is not None else None,
        "scl_row": sample.row,
        "scl_col": sample.col,
        "scl_asset_key": sample.asset_key,
        "scl_asset_href": sample.asset_href,
        "sampling_error": sample.error,
    }


def date_from_record(record: dict[str, Any]) -> str | None:
    for key in ("pixel_acquisition_time_utc", "granule_sensing_time_utc", "acquisition_time_utc"):
        parsed = parse_datetime(record.get(key))
        if parsed is not None:
            return parsed.date().isoformat()
    return None


def simple_csv_record(record: dict[str, Any]) -> dict[str, Any]:
    return {
        "product_id": record.get("item_id"),
        "latitude": record.get("latitude"),
        "longitude": record.get("longitude"),
        "date_utc": date_from_record(record),
        "scene_time_utc": record.get("acquisition_time_utc"),
        "granule_time_utc": record.get("granule_sensing_time_utc"),
        "pixel_time_utc": record.get("pixel_acquisition_time_utc"),
        "cloud_cover_percent": record.get("item_cloud_cover_percent"),
        "pixel_is_cloud": record.get("pixel_is_cloud"),
        "scl_value": record.get("scl_value"),
        "scl_class": record.get("scl_label"),
        "band_number": record.get("pixel_acquisition_time_band_id"),
    }


def write_csv(path: str, records: list[dict[str, Any]]) -> None:
    fieldnames = [
        "product_id",
        "latitude",
        "longitude",
        "date_utc",
        "scene_time_utc",
        "granule_time_utc",
        "pixel_time_utc",
        "cloud_cover_percent",
        "pixel_is_cloud",
        "scl_value",
        "scl_class",
        "band_number",
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(simple_csv_record(record) for record in records)


def write_jsonl(path: str, records: list[dict[str, Any]]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=False) + "\n")


def main() -> int:
    args = parse_args()
    if not -90 <= args.lat <= 90:
        raise SystemExit("--lat must be between -90 and 90.")
    if not -180 <= args.lon <= 180:
        raise SystemExit("--lon must be between -180 and 180.")

    records: list[dict[str, Any]] = []
    token = resolve_cdse_token(args)
    if needs_cdse_token(args) and not token:
        raise SystemExit(
            "No CDSE access token was provided. Set CDSE_ACCESS_TOKEN, set CDSE_USERNAME/CDSE_PASSWORD, "
            "pass --cdse-token, or rerun without --no-auth-prompt."
        )
    for item in iter_stac_items(args):
        sample = (
            PixelSample(None, None, None, None, None)
            if args.no_sample
            else sample_scl(item, args.lon, args.lat, token, args.cache_dir)
        )
        detector_sample = resolve_detector_id(item, args, token)
        timing = estimate_pixel_timing(
            item,
            sample,
            args.band_id,
            detector_sample.detector_id,
            args.metadata_root,
            token,
            args.timing_origin,
            args.time_per_line,
            args.detector_parallax,
            args.granule_sensing_time_position,
            detector_sample.source,
        )
        if detector_sample.error and not timing.error:
            timing.error = f"Detector sampling warning: {detector_sample.error}"
        elif detector_sample.error and timing.error:
            timing.error = f"{timing.error} Detector sampling warning: {detector_sample.error}"
        records.append(item_record(item, sample, timing, args.lon, args.lat))
        print(
            f"{len(records):5d} {records[-1]['acquisition_time_utc']} "
            f"{records[-1]['item_id']} pixel_time={records[-1]['pixel_acquisition_time_utc'] or records[-1]['pixel_acquisition_time_min_utc'] or 'n/a'} "
            f"SCL={records[-1]['scl_value']} "
            f"{records[-1]['scl_label'] or records[-1]['sampling_error'] or 'not sampled'}",
            file=sys.stderr,
        )

    write_csv(args.out, records)
    if args.jsonl:
        write_jsonl(args.jsonl, records)

    print(f"Wrote {len(records)} records to {args.out}")
    if args.jsonl:
        print(f"Wrote JSONL records to {args.jsonl}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
