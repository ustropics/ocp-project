# set  general variables
year = 2024
deg = 25

# get prefix based on degree given
if deg == 25:
    deg_prefix = 'GOMb0.04'

# set variables for hycom data retrieval 
dir_data_2d = f'output/data/{deg}deg/2d/'
dir_data_3d = f'output/data/{deg}deg/3d/'

dir_csv = f'output/csv/'



url_catalog = f'https://tds.hycom.org/thredds/catalog/datasets/{deg_prefix}/reanalysis/data/daily_netcdf/{year}/catalog.html'
url_data = f'https://tds.hycom.org/thredds/fileServer/datasets/{deg_prefix}/reanalysis/data/daily_netcdf/{year}/'