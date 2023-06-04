"""matching_clusters.py

Created on Fri May 26 15:10:11 2023

@author: McK
last edited: Sunday June 4, 2023

now matches clusters within 1 deg of each other without first using topcat's join method.

writes region files of each lovoccs cluster with its potential matches as the coordinates
of the region file. these matches are based off of calculating the angular distance
between the lovoccs cluster and the potential match. if the distance is less than
2 times the r500 of the lovoccs cluster, then i kept the match.
"""
#my imports
import pandas as pd
import numpy as np
from astropy import units as u
import csv
import subprocess
import os
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table

#imports from David Turner
from astropy.cosmology import Cosmology
from astropy.units import Quantity
from astropy.cosmology import LambdaCDM
DEFAULT_COSMO = LambdaCDM(70, 0.3, 0.7)

#importing data from big catalog and lovoccs csv file
desi = fits.open(r'/Users/McK/Desktop/galaxy_clusters_desidr9.fits')
lovoccs = pd.read_csv(r'/Users/McK/Desktop/comb_ims_for_mckenna/lovoccs_xmm_subsamp.csv')
table = Table(desi[1].data)

ra_desi = table['RA_PEAK']
dec_desi = table["DEC_PEAK"]
z_desi = table["PHOTO_Z_PEAK"]
r500_desi = table["R_500"]

ra_lovoccs = lovoccs["ra"]
dec_lovoccs = lovoccs["dec"]
z_lovoccs = lovoccs["redshift"]
r500_lovoccs = lovoccs["R500"]
name_lov = lovoccs["name"]

catalog1_coords = SkyCoord(ra=ra_lovoccs * u.deg, dec=dec_lovoccs * u.deg)
catalog2_coords = SkyCoord(ra=ra_desi * u.deg, dec=dec_desi * u.deg)

#empty list to store matches
matched_clusters = []

#loop through lovoccs clusters
for i, cluster1_coord in enumerate(catalog1_coords):
    separations = cluster1_coord.separation(catalog2_coords) #calculate angular separation with catalog2 clusters
    within_1deg = np.where(separations < 1 * u.deg)[0] #find catalog2 clusters within 1 degree

    if len(within_1deg) > 0:
        matching_indices = within_1deg.tolist() #create a list of matching cluster indices
        matched_clusters.append([i] + matching_indices) #append the list to the matched_clusters list
        
total_matches = sum(len(matches[1:]) for matches in matched_clusters)
#same number as before using TOPCAT

#getting redshift, r500, ra, dec vals from matches
z = [] #redshift for lovoccs cluster
r500 = [] #r500 of lovoccs cluster
ra_lov = [] #RA of lovoccs cluster
dec_lov = [] #dec of lovoccs cluster
ra_other = [] #RA from other cluster
dec_other = [] #dec from other cluster
r500_other = [] #r500 of matched cluster
z_other = [] #redshift of matched cluster
name = [] #name of lovoccs cluster

for matches in matched_clusters:
    catalog1_index = matches[0] #index of lovoccs cluster
    catalog2_indices = matches[1:] #index of catalog cluster that matches lovoccs
    
    for catalog2_index in catalog2_indices:
        
        #for each macthed cluster, append the following:
        ra_lov.append(ra_lovoccs[catalog1_index])
        ra_other.append(ra_desi[catalog2_index])
        dec_lov.append(dec_lovoccs[catalog1_index])
        dec_other.append(dec_desi[catalog2_index])
        z.append(z_lovoccs[catalog1_index])
        z_other.append(z_desi[catalog2_index])
        r500.append(r500_lovoccs[catalog1_index])
        r500_other.append(r500_desi[catalog2_index])
        name.append(name_lov[catalog1_index])

def rad_to_ang(rad: Quantity, z: float, cosmo: Cosmology = DEFAULT_COSMO) -> Quantity:
    """
    Written by David Turner
    Editted by McKenna Leichty
    
    Converts radius in length units to radius on sky in degrees.
    :param Quantity rad: Radius for conversion.
    :param Cosmology cosmo: An instance of an astropy cosmology, the default is a flat LambdaCDM concordance model.
    :param float z: The _redshift of the source.
    :return: The radius in degrees.
    :rtype: Quantity
    """
    distance_Mpc = rad*u.Mpc
    d_a = cosmo.angular_diameter_distance(z)
    ang_rad = (distance_Mpc / d_a) * (180 / np.pi)
    return ang_rad #in degrees

#1. r500 of lovoccs cluster converted to degrees
r500_deg = np.zeros(len(r500))

for i in range(len(r500_deg)):
    r500_deg[i] = rad_to_ang(r500[i], z[i])

#2. calculate the distance between lovoccs cluster and matched cluster
def ang_distance(ra_lov, ra_other, dec_lov, dec_other):
    """
    calculates the distance between lovoccs cluster and matched cluster in degrees.

    Parameters
    ----------
    ra_lov : RA of lovoccs cluster in decimal degrees
    ra_other : RA of matched cluster in decimal degrees
    dec_lov : dec of lovoccs cluster in decimal degrees
    dec_other : dec of matched cluster in decimal degrees

    Returns
    -------
    Angular distance between clusters in degrees.

    """
    ang_distance = np.sqrt(((ra_lov - ra_other)*np.cos(dec_lov))**2 + (dec_lov - dec_other)**2)
    return ang_distance #in decimal degrees

ang_distances = np.zeros(len(r500))
for i in range(len(ang_distances)):
    ang_distances[i] = ang_distance(ra_lov[i], ra_other[i], dec_lov[i], dec_other[i])

#3. find indices of real matches
def matches(r500_deg, ang_distances):
    """
    if distances are less than 2 times the r500 of lovoccs cluster, keep 
    index of these matches.

    Parameters
    ----------
    r500_deg : r500 of lovoccs cluster in degrees
    ang_distances : distance between lovoccs and matched cluster in degrees

    Returns
    -------
    idx : the index of matches in input.csv to keep

    """
    idx = []
    for i in range(len(r500_deg)):
        if ang_distances[i] < 2*r500_deg[i]:
            idx.append(i)
    return idx

kept_matches = matches(r500_deg, ang_distances)

#4. store these rows with matches in them into a new csv file (not really needed)
def write_file(input_file, output_file, rows_to_keep):
    """
    writes a new csv file that keeps the rows of real matches based on indices just found
    in kept_matches list.

    Parameters
    ----------
    input_file : name of input file with all matches
    output_file : name of file you want created with good matches
    rows_to_keep : indices of rows you want to keep in output file

    Returns
    -------
    output_file with all good matches

    """
    with open(input_file, 'r') as csv_input, open(output_file, 'w', newline='') as csv_output:
        reader = csv.reader(csv_input)
        writer = csv.writer(csv_output)
        for row_index, row in enumerate(reader):
            if row_index == 0 or (row_index - 1) in rows_to_keep:
                writer.writerow(row)
    csv_input.close()
    csv_output.close()

rows_to_keep = kept_matches
input_file = r'/Users/McK/Desktop/desidr9_matches.csv'
output_file = r'/Users/McK/Desktop/test_desidr9_matches.csv'
#write_file(input_file, output_file, rows_to_keep)

#5. separate each lovoccs cluster with its match
#we will end up with a list of lists. i.e. lov_match[0] will be matched clusters idx with first lovoccs cluster
lov_match = []
idx_ls = [] #start off with empty list
for i in range(len(kept_matches)):
    if i == len(kept_matches)-1: #if at end of list, add last idx to idx_ls
        idx_ls.append(kept_matches[i])
        lov_match.append(idx_ls) #all idx for one lovoccs cluster appended to big list
    
    elif ra_lov[kept_matches[i]] == ra_lov[kept_matches[i+1]]: #if same lovoccs cluster,
        idx_ls.append(kept_matches[i]) #append idx to ls
    
    else: #if not same lovoccs cluster, append to main list and start over
        lov_match.append(idx_ls) #all idx for one lovoccs cluster appended to big list
        idx_ls = [] #reset list

#92 matched lovoccs clusters: len(lov_match)

#6. create coordinates for each lovoccs cluster
def coordinate(lov_match):
    """
    writes coordinates for each lovoccs clusters' matched clusters in (ra, dec, radius) form.

    Parameters
    ----------
    lov_match : lists of each lovoccs clusters' matched clusters

    Returns
    -------
    coord : list of coordinates for each lovoccs cluster

    """
    coord = []
    for idx in lov_match:
        ra = ra_other[idx]
        dec = dec_other[idx]
        rad = rad_to_ang(r500_other[idx], z_other[idx]) #in degrees
        
        #convert to arcseconds
        rad_arcsec = rad * 3600 #1 deg * (60 arcmin / 1 deg) * (60 arcsec / 1 arcmin)
        cluster = (ra, dec, rad_arcsec)
        coord.append(cluster)
    
    return coord

lov_coords = []
for i in range(len(lov_match)): #for each lovoccs cluster, write coordinates
    lov_coord = coordinate(lov_match[i])
    lov_coords.append(lov_coord)
    
#7. writing region files for each lovoccs cluster and its matches
def region_file(coords, name, directory):
    """
    Writes a region file with coordinates of matched clusters.

    Parameters
    ----------
    coords : list of values related to matched cluster;
             ex: (ra, dec, radius")
    name : name of the lovoccs cluster

    Returns
    -------
    Saves a region file.

    """
    filename = os.path.join(directory, name + '.reg')  # Construct the full file path using the directory and cluster name
    file_content = '# Region file format: DS9 version 4.1\n'
    file_content += 'global color=red dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
    file_content += 'fk5\n'
    
    for coordinate in coords:
        if len(coordinate) == 3:
            file_content += f'circle({coordinate[0]}d,{coordinate[1]}d,{coordinate[2]}")\n'
        else:
            continue
    with open(filename, 'w') as file:
        file.write(file_content)

#name = matched_table["name"] #name of lovoccs cluster
dr = r'/Users/McK/Desktop/reu2023/ds9_test'
#already saved files to above directory
#for i in range(len(lov_coords)):
    #print(i)
    #if i == 31:
        #continue
    #else:
        #region_file(lov_coords[i], name[lov_match[i][0]], dr)

#8. opening ds9, uploading the fits and region files, saving image to output directory
def open_fits_with_region(fits_file, region_file, output_directory):
    """
    runs on the command line to activate environment where ds9 is installed,
    opens fits file on ds9, uploads region file, changes scale to log and 
    color to viridis, saves image as .png, then exits ds9.

    Parameters
    ----------
    fits_file : directory to fits file
    region_file : directory to region file
    output_directory : directory where images should be saved

    Returns
    -------
    Saves .png images to output_directory.

    """
    conda_env_name = "ciao-4.15"

    # Construct the command to open ds9, save the image, and quit ds9
    output_filename = os.path.basename(fits_file).replace('.fits', '_overlay.png')
    output_path = os.path.join(output_directory, output_filename)

    activate_env = f"conda run -n {conda_env_name} ds9 {fits_file} -regions load {region_file} -scale log -cmap viridis -saveimage {output_path} -quit"

    # Run the command
    subprocess.call(activate_env, shell=True)

def process_fits_files(directory, output_directory):
    """
    goes through each fits file in input directory, if there is a matching .reg
    file name then call on open_fits_with_region(). if there is not a matching
    region file, then print that there was no region file found and skip to next
    fits file.

    Parameters
    ----------
    directory : directory path where fits and region files can be found
    output_directory : directory path where images should be saved

    Returns
    -------
    saves .png images in output_directory

    """
    fits_files = [file for file in os.listdir(directory) if file.endswith('.fits')]

    for fits_file in fits_files:
        # Generate the corresponding region file name
        region_file = fits_file.replace('.fits', '.reg')

        # Check if the region file exists
        if os.path.exists(os.path.join(directory, region_file)):
            # Call the open_fits_with_region function with the output directory
            open_fits_with_region(os.path.join(directory, fits_file), os.path.join(directory, region_file), output_directory)
        else:
            print(f"No region file found for {fits_file}. Skipping...")

input_directory = r'/Users/McK/Desktop/reu2023/ds9_test'
output_directory = r'/Users/McK/Desktop/reu2023/ds9_test/saved_images'
#process_fits_files(input_directory, output_directory)

""" these are the fits files which had no region files (25 out of 116 files):
    
No region file found for MCXCJ0910.6-1034.fits. Skipping...
No region file found for MCXCJ1133.2+6622.fits. Skipping...
No region file found for MCXCJ1539.5-8335.fits. Skipping...
No region file found for MCXCJ1347.4-3250.fits. Skipping...
No region file found for MCXCJ2347.7-2808.fits. Skipping...
No region file found for MCXCJ1254.6-2913.fits. Skipping...
No region file found for MCXCJ1326.9-2710.fits. Skipping...
No region file found for MCXCJ2331.2-3630.fits. Skipping...
No region file found for MCXCJ1257.2-3022.fits. Skipping...
No region file found for MCXCJ1244.6-1159.fits. Skipping...
No region file found for MCXCJ0013.6-1930.fits. Skipping...
No region file found for MCXCJ2313.0-2137.fits. Skipping...
No region file found for MCXCJ1257.1-1724.fits. Skipping...
No region file found for MCXCJ2312.3-2130.fits. Skipping...
No region file found for MCXCJ1327.9-3130.fits. Skipping...
No region file found for MCXCJ0257.8+1302.fits. Skipping...
No region file found for MCXCJ0258.9+1334.fits. Skipping...
No region file found for MCXCJ1141.4-1216.fits. Skipping...
No region file found for MCXCJ0909.1-0939.fits. Skipping...
No region file found for MCXCJ2152.4-1933.fits. Skipping...
No region file found for MCXCJ2316.1-2027.fits. Skipping...
No region file found for MCXCJ2034.7-3548.fits. Skipping...
No region file found for MCXCJ1130.3-1434.fits. Skipping...
No region file found for MCXCJ1558.3-1410.fits. Skipping...
No region file found for MCXCJ2149.1-3041.fits. Skipping...
"""





