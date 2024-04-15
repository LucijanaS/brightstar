# brightstar
This repository uses the Yale Bright Star Catalog which contains 9110 of the brightest stars
(found at http://tdc-www.harvard.edu/catalogs/bsc5.html).
The cataloge is in ASCII format and was converted to a .JSON format in the repository 
https://github.com/brettonw/YaleBrightStarCatalog.
The data file used from that repository is bsc5-all.json.

## brightstar_input.py
Here one should input the locations(s) of observation, the date, the UTC offset, number of brightest stars one is interested in, and so on.

## brightstar_functions.py
All functions used in the following Python scripts can be found here.

## 1create_list.py
With 1create_list.py the .json file will be opened, scanned through for the properties important for our purposes, and the 
n brightest stars are saved in the file nstars_data.csv, where n is specified in brightstar_input.

## 2find_stars.py
Now instead of going through all the stars saved in 1000stars_data.csv of them to find the stars visible on a certain night, one can use the python script 
2find_stars.py.
It will use RA and Dec of all stars to determine if a star is visible for most of the night. Time of observation was assumed to be between 8pm and 8am of the day chosen, and a star counts as visible in that night if the altitude is >10° for at least 3/4 of that night. 
A list of visible stars is then be saved as stars_visible_date.csv file. 

## 3filter_stars.py
If the list of candidate stars is still too extensive one can filter out by magnitude, estimated apparent size and spectral photon flux density.
A new csv file will be created if parameters specified.

## 4UVW_plane.py
For the case where intensity inteferrometry is done at two different locations of known coordinates that were put in brightstar_input.py, one can determine the projection of the star on the UVW-plane during the night. Using the projection together with the squared visibility in that plane one can determine the stars which would cover the most area in the UVW-plane and the most intensity variation, leading to a more accurate measurement of the diameter of the star. The values will be added to the previous .csv file and the stars with most intensity variations will be plotted.

## 5bright_star.py
If a suitable star was found, one can use the RA and Dec of that star and the diameter or input the row from an interesting star extracted from the stars_visible_date.csv file and to plot the star's position on the night sky in horizontal coordinates as seen from the observers location. Further it plots the UVW-plane projection for that particular night, star and location together with the intensity pattern.

## main.py
Run main.py to run the first 4 python scripts.

## brightest_stars.ipynb
The Jupyter notebook brightest_stars.ipynb explains some of the details of how for example the diameter of the star is
determined based on its effective temperature and the apparent magnitudes in different wavelengths.

## SNR_test.ipynb
The Jupyter notebook SNR_test.ipynb is there to experiment with input values of efficiency, spectral photon flux density, electric bandwith from which one can try to estimate what diameter/area and observation time to gather significant results.
