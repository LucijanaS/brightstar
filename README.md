# brightstar
This repository uses the Yale Bright Star Catalog which contains 9110 of the brightest stars
(found at http://tdc-www.harvard.edu/catalogs/bsc5.html).
The cataloge is in ASCII format and was converted to a .JSON format in the repository 
https://github.com/brettonw/YaleBrightStarCatalog.
The data file use is bsc5-all.json.

## brightstar_input.py
Here one should input the locations(s) of observation, the date, the UTC offset

## brightstar_functions.py
All functions used in the following Python scripts can be found here.

## 1create_list.py
With 1create_list.py the .json file was opened, scanned through for the properties important for our purposes, and the 
1000 brightest stars were saved in the file 1000stars_data.csv. This code can be adapted to create list with more or
less stars.

## 2find_stars.py
Now instead of going through all the stars saved in 1000stars_data.csv of them to find the stars with suiting properties, one can use the python script 
2find_stars.py.
It will use RA and Dec of all stars to determine if a star is visible for most of the night (time of observation was assumed to be between 8pm and 8am of the day chosen). 
A list of visible stars can then be saved as a csv. file. If the list of candidate stars is still too extensive one can
look for stars with a certain diameter and create another list and .csv of those.

## 3UVW_plane.py
For the case where intensity inteferrometry is done at two different locations of known coordinates that were put in brightstar_input.py, one can determine the projection of the star on the UVW-plane during the night. Using the projection together with the intensity pattern in that plane one can determine the stars which would cover the most area in the UVW-plane and the most intensity variation, leading to a more accurate measurement of the diameter of the star. The values will be added to the previous .csv file and the stars with most intensity variations will be plotted.

## 4bright_star.py
If a suitable star was found, one can use the RA and Dec of that star and the diameter and the code 4bright_star.py to plot the star's 
position on the night sky in horizontal coordinates as seen from the observers location (again location and date need 
to be entered). Further it computes the UVW-plane projection for that particular night, star and location.

## main.py
Run main.py to run the first 3 python scripts.

## brightest_stars.ipynb
The Jupyter notebook brightest_stars.ipynb explains some of the details of how for example the diameter of the star is
determined based on its effective temperature and the apparent magnitudes in different wavelengths.
