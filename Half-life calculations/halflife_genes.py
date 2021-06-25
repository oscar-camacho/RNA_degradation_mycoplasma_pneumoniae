import numpy as np
import glob
import re
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os

def coverage_normalization():
    """Coverage of each CDS is normalized dividing by the sum of the coverages of the 3 rRNAs"""
    #Normalize coverage in each sample:
    filenames = glob.glob("*average_coverage.bed")
    file_dict = dict()
    for filename in filenames:
        bed_file = open(filename, "r")
        coverage_list = []
        common_information_list = []
        for line in bed_file:
            match = re.match(r"(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)", line)
            coverage = float(match.group(7))
            coverage_list.append(coverage)
            if match.group(4) == "Mpnr01":
                ribosomal_cov1 = float(match.group(7))
            if match.group(4) == "Mpnr02":
                ribosomal_cov2 = float(match.group(7))
            if match.group(4) == "Mpnr03":
                ribosomal_cov3 = float(match.group(7))
            if filename == filenames[-1]:
                common_information = match.group(1) + "\t" + match.group(2) + "\t" + match.group(3) + "\t" + match.group(4) + "\t" + match.group(5) + "\t" + match.group(6)
                common_information_list.append(common_information)

        sum_ribosomal_cov = ribosomal_cov1 + ribosomal_cov2 + ribosomal_cov3
        normalized_coverage_list = [k/sum_ribosomal_cov for k in coverage_list]
        file_dict[filename] = normalized_coverage_list

    ## Write merging results:
    merging_file = open("merged_coverage.bed", "w")
    # Sorting
    def atoi(text):
        return int(text) if text.isdigit() else text

    def natural_keys(text):
        return [ atoi(c) for c in re.split(r'(\d+)', text) ]
    filenames.sort(key=natural_keys)

    # Writing header
    merging_file.write("\t\t\t\t\t\t")
    for filename in filenames:
        match = re.match(r"^(Nub_\d+.?)_", filename)
        mini_name = match.group(1)
        merging_file.write(mini_name)
        if filename != filenames[-1]:
            merging_file.write("\t")
        else:
            merging_file.write("\n")

    for common_information in common_information_list:
        merging_file.write(common_information + "\t")
        for file in filenames:
            merging_file.write(str(file_dict[file][0]))
            del file_dict[file][0]
            if file != filenames[-1]:
                merging_file.write("\t")
            else:
                merging_file.write("\n")

coverage_normalization()


def fitting_regression():
    file = open("merged_coverage.bed", "r")
    new_file = open("halflife_genes.bed", "w")
    new_file.write("Transcript\tinitial_coord\tfinal_coord\tstrand\tk(decay)\thalf_live\tr-quared\n")
    x_data = []
    filenames = glob.glob("*average_coverage.bed")
    for filename in filenames:
        match = re.match(r"^Nub_(\d+)._*", filename)
        number = int(match.group(1))
        x_data.append(number)
    x_data.sort()

    x_data = np.array(x_data)
    y_data = []
    for line in file:
        if line.startswith("\t"):
            continue
        else:
            elements = line.split()
            transcript = elements[3]
            initial_coord = elements[1]
            final_coord = elements[2]
            y_list = elements[6:]
            y_list = list(map(float, y_list))
            for number in y_list:
                if number == 0:
                    index = y_list.index(number)
                    y_list[index] = 0.0000002

            y_data = np.array(y_list)
            log_y_data = np.log(y_data)
            curve_fit = np.polyfit(x_data, log_y_data, 1)
            k_decay, c = curve_fit
            p = np.poly1d(curve_fit)
            yhat = p(x_data)
            ybar = sum(log_y_data)/len(log_y_data)
            SST = sum((log_y_data - ybar)**2)
            SSreg = sum((yhat - ybar)**2)
            R2 = SSreg/SST
            #print("R-squared = " + str(R2))
            # Calculate half-lives:
            half_live = np.log(2) / (-k_decay)
            #print("Transcript: " + transcript + "\t" + initial_coord + "-" + final_coord + "\nK decay=" + str(-k_decay) + "\tc=" + str(c) + "\tHalf-live=" + str(half_live))

            # Visalizing data:
            y = np.exp(k_decay*x_data) * np.exp(c)
            plt.plot(x_data, y_data, "o")
            plt.plot(x_data, y)
            #plt.show()

            new_file.write(transcript + "\t" + initial_coord + "\t" + final_coord + "\t" + elements[5] + "\t" + str(-k_decay) + "\t" + "\t" + str(half_live) + "\t" + str(R2) + "\n")




fitting_regression()
