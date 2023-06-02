# Background_Clusters

From this website, (https://www.scidb.cn/en/detail?dataSetId=7797b553a23846a187b7746d8fc555a5) I downloaded three files (galaxy_clusters_desdr2.fits, galaxy_clusters_desidr9.fits, and galaxy_clusters_hscpdr3_wide.fits) to match the nearest clusters from these data to the LOVOCCS clusters. 

To start, I used the program TOPCAT to match the big catalogue clusters with the LOVOCCS clusters within 1 degree on the sky. These new matches were saved as a .csv file from which I used in my matching_clusters.py code. In this Python script, I converted the LOVOCCS clusters r500 to decimal degrees and calculated the distance (in decimal degrees) from the LOVOCCS clusters and their matches. If these distances were less than 2 times the r500 of the LOVOCCS cluster then I considered these a "true" match and saved them. These "true" matches can be written into a new .csv file using the write_file() function. 

From here, I wrote .reg region files for each LOVOCCS cluster and its "true" matches. First, I separated each LOVOCCS cluster and its matches into a separate list, and then created coordinates for the matched clusters per LOVOCCS cluster. The coordinates were stored in the form: circle(RAd, decd, radius") where the RA and dec are in decimal degrees and the r500 of the matched cluster was stored in arcseconds.

Using these coordinates and the function region_file(), I wrote .reg files for each LOVOCCS cluster and its matches in the fk5, ds9 format. These region files were saved in a specified directory.

Finally, I used the subprocess.call() function to call the command line from inside the Python script. I had a specified environment where my ds9 command was installed, so I first activated the conda environment. The ds9 was called to open the fits file, read in the corresponding .reg file (based on if the .reg file had the same name as the fits file), changed the scale to log and the color to viridis, saved the image as a .png to a specified directory, and then closed ds9.

This was repeated for every fits image in the specified directory. I had some fits files without a .reg file so if there was no corresponding .reg file, the fits file was skipped and printed out the fits file name.
