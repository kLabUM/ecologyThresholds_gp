# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "mcp[cli]>=1.2.0",
#   "httpx>=0.27",
# ]
# ///
"""
Minimal InfluxDB Cloud (v2 API) MCP server.

Exposes four tools:
  - list_buckets()           : list all buckets in the configured org
  - schema_explore(bucket)   : measurements, tag keys, and field keys for a bucket
  - query_flux(query)        : run any Flux query, return raw CSV
  - query_json(query)        : run any Flux query, return JSON records (one per row)

Reads three env vars (set in the Claude Desktop mcpServers config):
  INFLUX_URL   e.g. https://us-east-1-1.aws.cloud2.influxdata.com
  INFLUX_TOKEN your read-only API token
  INFLUX_ORG   e.g. "Digital Water Lab"
"""

import csv
import io
import json
import os
from typing import Any

import httpx
from mcp.server.fastmcp import FastMCP

INFLUX_URL = os.environ.get("INFLUX_URL", "").rstrip("/")
INFLUX_TOKEN = os.environ.get("INFLUX_TOKEN", "")
INFLUX_ORG = os.environ.get("INFLUX_ORG", "")

if not (INFLUX_URL and INFLUX_TOKEN and INFLUX_ORG):
    raise RuntimeError(
        "Missing one of INFLUX_URL / INFLUX_TOKEN / INFLUX_ORG environment variables."
    )

mcp = FastMCP("influxdb")


def _auth_headers(extra: dict[str, str] | None = None) -> dict[str, str]:
    headers = {"Authorization": f"Token {INFLUX_TOKEN}"}
    if extra:
        headers.update(extra)
    return headers


def _run_flux(query: str) -> str:
    """POST a Flux query to /api/v2/query, return CSV text. Raises on HTTP error."""
    r = httpx.post(
        f"{INFLUX_URL}/api/v2/query",
        params={"org": INFLUX_ORG},
        headers=_auth_headers({
            "Content-Type": "application/vnd.flux",
            "Accept": "application/csv",
        }),
        content=query,
        timeout=60.0,
    )
    if r.status_code != 200:
        raise RuntimeError(f"Influx query failed ({r.status_code}): {r.text[:500]}")
    return r.text


def _csv_to_records(csv_text: str) -> list[dict[str, Any]]:
    """
    Convert Influx annotated CSV into a flat list of records.

    Influx returns multi-table CSV with annotation rows (#datatype, #group, #default)
    and blank lines between tables. We strip annotations and concatenate all data rows
    from all tables, keyed by their header row.
    """
    records: list[dict[str, Any]] = []
    # Split on blank lines (table separators) then parse each table independently
    tables = [t for t in csv_text.split("\r\n\r\n") if t.strip()]
    if len(tables) == 1:
        tables = [t for t in csv_text.split("\n\n") if t.strip()]
    for table in tables:
        lines = [ln for ln in table.splitlines() if ln.strip()]
        # Drop annotation rows (start with '#')
        data_lines = [ln for ln in lines if not ln.startswith("#")]
        if len(data_lines) < 2:
            continue
        reader = csv.DictReader(io.StringIO("\n".join(data_lines)))
        for row in reader:
            # Drop empty leading column that Influx CSV includes
            row.pop("", None)
            records.append(row)
    return records


@mcp.tool()
def list_buckets() -> str:
    """List all buckets in the configured InfluxDB org. Returns JSON."""
    r = httpx.get(
        f"{INFLUX_URL}/api/v2/buckets",
        params={"org": INFLUX_ORG, "limit": 100},
        headers=_auth_headers(),
        timeout=30.0,
    )
    if r.status_code != 200:
        raise RuntimeError(f"list_buckets failed ({r.status_code}): {r.text[:500]}")
    data = r.json()
    # Slim down to the useful fields
    buckets = [
        {"id": b.get("id"), "name": b.get("name"), "retentionRules": b.get("retentionRules")}
        for b in data.get("buckets", [])
    ]
    return json.dumps(buckets, indent=2)


@mcp.tool()
def schema_explore(bucket: str) -> str:
    """Return measurements, tag keys, and field keys for the given bucket. Returns JSON."""
    result: dict[str, list[str]] = {}
    queries = {
        "measurements": (
            'import "influxdata/influxdb/schema"\n'
            f'schema.measurements(bucket: "{bucket}")'
        ),
        "tagKeys": (
            'import "influxdata/influxdb/schema"\n'
            f'schema.tagKeys(bucket: "{bucket}")'
        ),
        "fieldKeys": (
            'import "influxdata/influxdb/schema"\n'
            f'schema.fieldKeys(bucket: "{bucket}")'
        ),
    }
    for name, q in queries.items():
        try:
            csv_text = _run_flux(q)
            records = _csv_to_records(csv_text)
            # schema.* queries return a "_value" column with each name
            result[name] = sorted({r.get("_value", "") for r in records if r.get("_value")})
        except Exception as e:
            result[name] = [f"ERROR: {e}"]
    return json.dumps(result, indent=2)


@mcp.tool()
def query_flux(query: str) -> str:
    """Run an arbitrary Flux query against the org. Returns the raw CSV response."""
    return _run_flux(query)


@mcp.tool()
def query_json(query: str) -> str:
    """Run an arbitrary Flux query against the org. Returns a JSON array of records."""
    csv_text = _run_flux(query)
    records = _csv_to_records(csv_text)
    return json.dumps(records)


if __name__ == "__main__":
    mcp.run()
