import argparse
import csv
import os

import requests

# ---------------- CONFIG ----------------
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SCHEME_DIR = os.path.join(BASE_DIR, "schemes")
BASE_URL = "https://lpc.web.cern.ch/fillingSchemes/2018"

# Load run_to_fill mapping
run_to_fill = {}
with open(os.path.join(BASE_DIR, "mono_run_to_fill.csv"), "r") as f:
  reader = csv.reader(f)
  for row in reader:
    try:
      run, fill = row[0].split(";")
      run = int(run.strip())
      fill = int(fill.strip())
      run_to_fill[run] = fill
    except ValueError:
      continue

# Load fill_to_scheme mapping
fill_to_scheme = {}
with open(os.path.join(BASE_DIR, "mono_fill_to_scheme.csv"), "r") as f:
  reader = csv.reader(f)
  for row in reader:
    try:
      fill, scheme = row[0].split(";")
      fill = int(fill.strip())
      scheme = scheme.strip()
      fill_to_scheme[fill] = scheme
    except ValueError:
      continue


# ------------- DOWNLOAD -------------
def ensure_scheme_file(scheme_name):
  os.makedirs(SCHEME_DIR, exist_ok=True)
  path = os.path.join(SCHEME_DIR, f"{scheme_name}.csv")

  if not os.path.exists(path):
    url = f"{BASE_URL}/{scheme_name}.csv"
    print(f"Downloading {url}")
    r = requests.get(url)
    r.raise_for_status()
    with open(path, "wb") as f:
      f.write(r.content)

  return path


# ------------- PARSER -------------
def load_scheme(scheme_name):
  path = ensure_scheme_file(scheme_name)

  beam1 = set()
  beam2 = set()

  with open(path, "r") as f:
    lines = f.readlines()

  # skip header line
  for line in lines[1:]:
    parts = [p.strip() for p in line.split(",")]

    if len(parts) < 4:
      continue

    try:
      ring = parts[2]
      rf_bucket = int(parts[3])
    except:
      continue

    # The LPC filling-scheme files store RF buckets (2.5 ns spacing), while the
    # rest of this script works with BX numbers (25 ns spacing, 0..3563).
    # Valid collision buckets are numbered 1, 11, 21, ... so we map them back
    # to BX indices with (rf_bucket - 1) / 10.
    bx = (rf_bucket - 1) // 10

    if ring == "ring_1":
      beam1.add(bx)
    elif ring == "ring_2":
      beam2.add(bx)

  colliding = beam1 & beam2

  return colliding


# ------------- CORE LOGIC -------------
_scheme_cache = {}


def get_colliding_for_run(run):
  fill = run_to_fill[run]
  scheme = fill_to_scheme[fill]

  if scheme not in _scheme_cache:
    _scheme_cache[scheme] = load_scheme(scheme)

  return _scheme_cache[scheme]


def has_prev_collision(run, bx):
  colliding = get_colliding_for_run(run)

  bx = bx % 3564
  prev = (bx - 1) % 3564

  return prev in colliding


def distance_to_prev_collision(run, bx, max_d=20):
  colliding = get_colliding_for_run(run)

  bx = bx % 3564

  for d in range(1, max_d + 1):
    if (bx - d) % 3564 in colliding:
      return d

  return None


def main():
  parser = argparse.ArgumentParser(
    description="Return distance to the previous collision for a given run and BX."
  )
  parser.add_argument("run", type=int, help="Run number")
  parser.add_argument("bx", type=int, help="BX number")
  parser.add_argument("max_distance", type=int, help="Maximum BX distance to search")
  args = parser.parse_args()

  distance = distance_to_prev_collision(args.run, args.bx, args.max_distance)
  print(distance if distance is not None else 0)


if __name__ == "__main__":
  main()
