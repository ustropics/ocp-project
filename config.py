# config.py (updated)

# Set general variables
year = 2024
day = '002'
deg = 25

# Get prefix based on degree given
if deg == 25:
    deg_prefix = 'GOMb0.04'

# Set variables for hycom data retrieval
dir_data_2d = f'output/data/{deg}deg/2d/'
dir_data_3d = f'output/data/{deg}deg/3d/'

dir_csv = f'output/csv/'

url_catalog = f'https://tds.hycom.org/thredds/catalog/datasets/{deg_prefix}/reanalysis/data/daily_netcdf/{year}/catalog.html'
url_data = f'https://tds.hycom.org/thredds/fileServer/datasets/{deg_prefix}/reanalysis/data/daily_netcdf/{year}/'
