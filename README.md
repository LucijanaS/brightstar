# brightstar
This repository uses the Yale Bright Star Catalog which contains 9110 of the brightest stars
(found at http://tdc-www.harvard.edu/catalogs/bsc5.html).
The cataloge is in ASCII format and was converted to a .JSON format in the repository 
https://github.com/brettonw/YaleBrightStarCatalog.
The data file use is bsc5-all.json.

With create_list.py the .json file was opened, scanned through for the properties important for our purposes, and the 
1000 brightest stars were saved in the file 1000stars_data.csv. This code can be adapted to create list with more or
less stars.

Now instead of going through all of them to find the stars with suiting properties, one can use the python script 
find_stars.py.
There one can enter the date of observation and location, and it will use RA and Dec of all stars to determine if a star
is visible for most of the night (time of observation was assumed to be between 8pm and 8am of the following day). 
A list of visible stars can then be saved as a csv. file. If the list of candidate stars is still too extensive one can
look for stars with a certain diameter and create another list and .csv of those.

If a suitable star was found, one can use the RA and Dec of that star and the code bright_star.py to plot the star's 
position on the night sky in horizontal coordinates as seen from the observers location (again location and date need 
to be entered). Further it computes the UVW-plane projection for that particular night, star and location.

The Jupyter notebook brightest_stars.ipynb explains some of the details of how for example the diameter of the star is
determined based on its effective temperature and the apparent magnitudes in different wavelengths.
