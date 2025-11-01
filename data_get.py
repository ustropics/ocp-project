from imports import *
from config import *


def get_3d_data():
    # Optional: Directory to save files (create if not exists)
    os.makedirs(dir_data_3d, exist_ok=True)

    # Fetch the catalog page
    response = requests.get(url_catalog)
    if response.status_code != 200:
        print(f"Failed to fetch catalog: {response.status_code}")
        exit(1)

    # Parse the HTML
    soup = BeautifulSoup(response.text, 'html.parser')

    # Find all <a> tags and extract file names ending with _3z.nc
    links = soup.find_all('a')
    file_names = [link.text.strip() for link in links if link.text.strip().endswith('_3z.nc')]

    if not file_names:
        print("No _3z.nc files found in the catalog.")
        exit(0)

    # Download each file if it doesn't exist locally
    for file_name in file_names:
        local_path = os.path.join(dir_data_3d, file_name)
        if os.path.exists(local_path):
            print(f"{file_name} already exists. Skipping download.")
            continue
        
        download_url = url_data + file_name
        print(f"Downloading {file_name} from {download_url}")
        
        # Stream the download for large files
        try:
            with requests.get(download_url, stream=True) as r:
                r.raise_for_status()
                with open(local_path, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
            print(f"Downloaded {file_name} successfully.")
        except requests.exceptions.RequestException as e:
            print(f"Failed to download {file_name}: {e}")

def get_2d_data():
    # Optional: Directory to save files (create if not exists)
    os.makedirs(dir_data_2d, exist_ok=True)

    # Fetch the catalog page
    response = requests.get(url_catalog)
    if response.status_code != 200:
        print(f"Failed to fetch catalog: {response.status_code}")
        exit(1)

    # Parse the HTML
    soup = BeautifulSoup(response.text, 'html.parser')

    # Find all <a> tags and extract file names ending with _3z.nc
    links = soup.find_all('a')
    file_names = [link.text.strip() for link in links if link.text.strip().endswith('_2d.nc')]

    if not file_names:
        print("No _2d.nc files found in the catalog.")
        exit(0)

    # Download each file if it doesn't exist locally
    for file_name in file_names:
        local_path = os.path.join(dir_data_2d, file_name)
        if os.path.exists(local_path):
            print(f"{file_name} already exists. Skipping download.")
            continue
        
        download_url = url_data + file_name
        print(f"Downloading {file_name} from {download_url}")
        
        # Stream the download for large files
        try:
            with requests.get(download_url, stream=True) as r:
                r.raise_for_status()
                with open(local_path, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
            print(f"Downloaded {file_name} successfully.")
        except requests.exceptions.RequestException as e:
            print(f"Failed to download {file_name}: {e}")
