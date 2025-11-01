from imports import *
from config import *

def eval_3d_data():

    if deg == 25:
        prefix = 'gomb4'
        
    # Hardcoded file path
    file_path = dir_data_3d + f'{prefix}_daily_2024_001_3z.nc'

    os.makedirs(dir_csv, exist_ok=True)
    csv_file = os.path.join(dir_csv, f"variables_{deg}deg_3d.csv")

    # Open the NetCDF file
    try:
        ds = nc.Dataset(file_path, 'r')
        print(f"Successfully opened: {file_path}\n")

        # Prepare data for CSV
        rows = []
        for var_name, var in ds.variables.items():
            dims = ', '.join(var.dimensions) if var.dimensions else ''
            shape = str(var.shape)
            dtype = str(var.dtype)

            # Get long_name and units from attributes
            long_name = var.getncattr('long_name') if 'long_name' in var.ncattrs() else ''
            units = var.getncattr('units') if 'units' in var.ncattrs() else ''

            rows.append({
                'variable': var_name,
                'dimensions': dims,
                'shape': shape,
                'data_type': dtype,
                'long_name': long_name,
                'units': units
            })

        # Write to CSV
        with open(csv_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=['variable', 'dimensions', 'shape', 'data_type', 'long_name', 'units'])
            writer.writeheader()
            writer.writerows(rows)

        print(f"Metadata saved to: {csv_file}")
        print(f"   Total variables exported: {len(rows)}")

        ds.close()

    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        print("   Make sure the file was downloaded using the previous script.")
    except Exception as e:
        print(f"Error reading file: {e}")


def eval_2d_data():

    if deg == 25:
        prefix = 'gomb4'
        
    # Hardcoded file path
    file_path = dir_data_2d + f'{prefix}_daily_2024_001_2d.nc'

    os.makedirs(dir_csv, exist_ok=True)
    csv_file = os.path.join(dir_csv, f"variables_{deg}deg_2d.csv")

    # Open the NetCDF file
    try:
        ds = nc.Dataset(file_path, 'r')
        print(f"Successfully opened: {file_path}\n")

        # Prepare data for CSV
        rows = []
        for var_name, var in ds.variables.items():
            dims = ', '.join(var.dimensions) if var.dimensions else ''
            shape = str(var.shape)
            dtype = str(var.dtype)

            # Get long_name and units from attributes
            long_name = var.getncattr('long_name') if 'long_name' in var.ncattrs() else ''
            units = var.getncattr('units') if 'units' in var.ncattrs() else ''

            rows.append({
                'variable': var_name,
                'dimensions': dims,
                'shape': shape,
                'data_type': dtype,
                'long_name': long_name,
                'units': units
            })

        # Write to CSV
        with open(csv_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=['variable', 'dimensions', 'shape', 'data_type', 'long_name', 'units'])
            writer.writeheader()
            writer.writerows(rows)

        print(f"Metadata saved to: {csv_file}")
        print(f"   Total variables exported: {len(rows)}")

        ds.close()

    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        print("   Make sure the file was downloaded using the previous script.")
    except Exception as e:
        print(f"Error reading file: {e}")