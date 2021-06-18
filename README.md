# SHARC - Version 2.1.06
#
      High Altitude Platform as IMT Base Station (HIBS) module.
	  - Author: Luciano Camilo Alexandre
	  - contact: luciano-camilo@hotmail.com

<p align="center">
  <img src="https://github.com/lucianodtel/sharc-TMGTelecom/blob/master/SHARC-development/sharc/img/sharc-logo.png?raw=true" width="550" title="SHARing and Compatibility studies of radiocommunication systems">
</p>

# Release notes:

#
	  - date: 18/06/2021
 - Refactoring ARNS.
#
	  - date: 17/06/2021
 - Refactoring ARNS.
 
#
	  - date: 31/05/2021

 - Added in /support the satellite/pseudo-satellite antenna beamwidth earth coverage script.
 - Added in /support the satellite/pseudo-satellite horizon visibility script.
 - Added in /support the antenna efficiency script.
 - HIBS Cluster options: 0 and 1;
 - HIBS Band 1 option;
 - Added 3D Antenna patterns for IMT Beamforming;
 - Added 3D Antenna patterns for Radars;
 - Added new Phased Array Radar antenna for ARNS and Meteorological Radar;
 - Added Pencil Radar antenna pattern;
 - Added Cosine^2 Radar antenna pattern;
 - Added 3D Pencil Radar antenna Plot;
 - Added 3D Cosine^2 Radar antenna Plot;
 - Added 3D Cosecant Squared Radar antenna plot;
 - Radar Phased Array Antenna element Cosine^2;
 - Added new section in parameters.ini regarding IMT Base Station;
 - Added new classes parameters_imt_base_station and parameters_antenna_imt_base_stations;
 - Added cossecant antenna pattern for Aeronautical Surveillance Radars;
 - Added cosine antenna pattern for Aeronautical Surveillance Radars according ITU-R M.1851-1 Table 2/3 and Fig. 4;
 - Added cossecant squared antenna pattern for Aeronautical Surveillance Radars according ITU-R M.1851-1 Figure 12;
 - Implementation of Earth-space path according ITU-R P.619 - Annex I-A;
 - Added new system: ARNS (Aeronautical Surveillance Radars and Metereological Radars);
 - Added new section in parameters.ini regarding ARNS station;
 - Added new system: SS-MLEO (Space station - Medium/Low earth orbit);
 - Added new section in parameters.ini regarding SS-MLEO station;
 - Added antenna pattern based on ITU-R S.1528-LEO (Annex 1);
 - Added azimuth and elevation statistical distributions in ARNS Section;

#
	  - date: 01/03/2021

 - Implemented Beamforming enabling/disabling option in parameters.ini;
 - PEP-8 Python code refactoring;
 - Added Bessel antenna pattern for Fixed Services;
 - Added ITU-R F.1245 Antenna Pattern for Fixed Services;
 - New antenna pattern based on ITU-R S.1528-LEO (Annex 1);
 - Added azimuth and elevation statistical distributions in Fixed Services Section;

#
	  - date: 20/12/2020

 - Power control (downlink) method for HIBS Base Station;
 - Added HIBS section in parameters.ini;
 - Added azimuth/elevation HIBS cell configuration.
 - Added conducted power per HIBS cell configuration.
 - Added back-off power per HIBS cell configuration.
 - Added different type of arrays in HIBS Base Station sectors;
 - Graphical module to import SHARC results in one plot;
 - Added HIBS Topology;
 - Added ITU-R F.1336 Antenna Pattern for Fixed Services;
 - Added azimuth and elevation statistical distributions in Radio Astronomy service section;
 - HIBS Topology (1, 3, 7 and 19 sectors);
 - Added in view.py HIBS section;
 - Added conducted power differences in HIBS Base Stations;

