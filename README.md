# TfL Crowding Data Download

Automated download of crowding profiles for **London Underground stations** from the [Transport for London Unified API](https://api.tfl.gov.uk/).

---

## What the data contains

The TfL Crowding API provides a **typical busyness profile** for each station, broken into 15-minute time bands across every day of the week. Values are expressed as a fraction of the station's historical peak (since July 2019), where `1.0` = the busiest the station has ever been.

This is a **relative** measure, not an absolute passenger count. It implicitly accounts for station size because the baseline is calibrated per station.

---

## Scripts

| Script | Description                                                                                          |
|---|------------------------------------------------------------------------------------------------------|
| `main.py` | Downloads 15-min crowding profiles for all LU stations, all 7 days. Outputs `tfl_crowding_data.csv`. |
| `summarise.py` | Aggregates the raw data into daily, weekly, and peak-hour summaries for further analysis.            |

---

## Setup

```bash
pip install requests pandas
```

An API key is **optional** but recommended to avoid rate limiting (~50 requests/min without a key). Get a free key at [api-portal.tfl.gov.uk](https://api-portal.tfl.gov.uk/).

Paste your key into `main.py`:

```python
API_KEY = "your_key_here"
```

---

## Usage

**1. Download the data**

```bash
python main.py
```

This may take some minutes. It produces `tfl_crowding_data.csv` (~170,000+ rows).

**2. Summarise**

```bash
python summarise.py
```

Produces three files:

| Output | Description |
|---|---|
| `tfl_crowding_daily.csv` | Mean / peak / min crowding per station per day of week |
| `tfl_crowding_weekly_avg.csv` | Single crowding value per station (averaged across all days) |
| `tfl_crowding_peak_hours.csv` | AM peak (07-09) and PM peak (17-19) crowding per station |

All outputs include `lat`, `lon`, and `naptan_id` for spatial joining.

---

## Output columns

### `tfl_crowding_data.csv`

| Column | Description |
|---|---|
| `station_name` | Station name |
| `naptan_id` | NaPTAN identifier (e.g. `940GZZLUOXC`) |
| `lat` | Latitude |
| `lon` | Longitude |
| `day_of_week` | Day of week (`Mon`, `Tue`, ..., `Sun`) |
| `time_band` | 15-min interval (e.g. `0800-0815`) |
| `hour` | Hour extracted from time band (0-23) |
| `pct_of_baseline` | Fraction of historical peak busyness (0.0 to 1.0) |

---

## Data source

- **API**: [TfL Unified API - Crowding endpoint](https://api.tfl.gov.uk/)
- **Licence**: [TfL Open Data](https://tfl.gov.uk/corporate/terms-and-conditions/transport-data-service) (Open Government Licence)
- **Baseline period**: July 2019 onwards
