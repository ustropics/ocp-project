from imports import *
from config import *
from concurrent.futures import ThreadPoolExecutor, as_completed
from bs4 import BeautifulSoup
import os
import requests


def download_file(file_name, dir_path, url_base):
    """Download a single file with streaming and error handling."""
    local_path = os.path.join(dir_path, file_name)
    if os.path.exists(local_path):
        print(f"{file_name} already exists. Skipping.")
        return file_name, True

    download_url = url_base + file_name
    print(f"Downloading {file_name} from {download_url}")

    try:
        with requests.get(download_url, stream=True, timeout=30) as r:
            r.raise_for_status()
            with open(local_path, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
        print(f"Downloaded {file_name} successfully.")
        return file_name, True
    except requests.exceptions.RequestException as e:
        print(f"Failed to download {file_name}: {e}")
        return file_name, False


def fetch_file_list(suffix, url_catalog):
    """Fetch and parse catalog, return list of filenames ending with suffix."""
    response = requests.get(url_catalog, timeout=30)
    if response.status_code != 200:
        raise RuntimeError(f"Failed to fetch catalog: {response.status_code}")

    soup = BeautifulSoup(response.text, 'html.parser')
    links = soup.find_all('a')
    return [link.text.strip() for link in links if link.text.strip().endswith(suffix)]


def download_dataset(suffix, dir_path, url_base, max_workers=6):
    """Download all files with given suffix in parallel."""
    os.makedirs(dir_path, exist_ok=True)

    try:
        file_names = fetch_file_list(suffix, url_catalog)
    except Exception as e:
        print(e)
        return

    if not file_names:
        print(f"No files ending with '{suffix}' found in catalog.")
        return

    print(
        f"Found {len(file_names)} files to download (max {max_workers} concurrently).")

    # Use ThreadPoolExecutor for parallel downloads
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_file = {
            executor.submit(download_file, fname, dir_path, url_base): fname
            for fname in file_names
        }

        for future in as_completed(future_to_file):
            fname = future_to_file[future]
            try:
                file_name, success = future.result()
                # Optional: track failed downloads
            except Exception as exc:
                print(f"{fname} generated an exception: {exc}")


# === Main functions ===
def get_3d_data():
    download_dataset('_3z.nc', dir_data_3d, url_data, max_workers=6)


def get_2d_data():
    download_dataset('_2d.nc', dir_data_2d, url_data, max_workers=6)
