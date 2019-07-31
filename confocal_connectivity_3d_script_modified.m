%July 30th, 2019
%This program analyzes the connectivity of water 
%Created by Geena.
%Modified by Allen and Edna.
%Note: this program makes use of the function "confocal_connectivity_3d"

%To use: from avizo export the excel files with data from the 'auto
%skeleton'and 'spatial graph statistics' modules for the lutidine and water 
%phases. Update lines 15-26 below and run the program. Output will be 
%formatted and saved in the file specified in "output_file".

clc
clear
warning off

%-->
series='Series026';
concentration='5050';
sample_name = ['Bipjels',series];
%-->
particles_skeleton_file = ['C:\Users\Ben\Documents\BIPJELS-Edna\5050\Test2and3\Series026\Binary\Bipjel.Smt.SptGraph.xml'];%From auto skeleton module
particles_statistics_file = ['C:\Users\Ben\Documents\BIPJELS-Edna\5050\Test2and3\Series026\Binary\Bipjel.Smt.statistics.xml'];%From spatial graph statistics module
output_file = ['Bipjels',series,'_Connectivity.xlsx'];

%Call function to calculate connectivity data. Return variable is a column
%vector
particles_data = confocal_connectivity_3d(particles_skeleton_file, particles_statistics_file);

%Write results to output file
headers = [{'Number of Nodes'};{'Max Coordination Number'}; ...
    {'Avg Coordination Number'};{'Standard Deviation Coordination Number'};{'Max Branch Length'}; ...
    {'Avg Branch Length'};{'Standard Deviation Branch Length'};{'Max Branch Diameter'};{'Avg Branch Diameter'};{'Standard Deviation Branch Diameter'}; ...
    {'Tails not on Boarder'};{'Tail Percentage'}];

%-->
cd 'C:\Users\Ben\Documents\BIPJELS-Edna\'

xlswrite(output_file1,[{'Clay'}],1,'B1:B1');
xlswrite(output_file1,headers,1,'A3:A12');
xlswrite(output_file1,[particles_data1],1,'B3:B12');

xlswrite(output_file2,[{'Water'}],1,'B1:B1');
xlswrite(output_file2,headers,1,'A3:A12');
xlswrite(output_file2,[particles_data2],1,'B3:B12');

%-->
xlswrite('Results_concentration.xlsx',particles_data',concentration,'I11:T11');

