import numpy as np
import glob
import re
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os



def calculate_normalization():
    pileup_filenames = glob.glob("*.strand*_coverage.bed")
    # Delete previous files:
    existing_files = glob.glob("*normalized_coverage.bed")
    for existing_file in existing_files:
        os.remove(existing_file)

    for pileup_filename in pileup_filenames:
        # Get normalization value using rRNAs in respective 'CDS_average_coverage.bed file'
        match = re.match(r"^(.*)strand._coverage.bed", pileup_filename)
        minifile_name = match.group(1)
        CDS_bed_file = open(minifile_name+"CDS_average_coverage.bed")
        for line in CDS_bed_file:
            line = line.strip()
            match = re.match(r"(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)", line)
            if match.group(4) == "Mpnr01":
                ribosomal_cov1 = float(match.group(7))
            if match.group(4) == "Mpnr02":
                ribosomal_cov2 = float(match.group(7))
            if match.group(4) == "Mpnr03":
                ribosomal_cov3 = float(match.group(7))
        sum_ribosomal_cov = ribosomal_cov1 + ribosomal_cov2 + ribosomal_cov3


        list = []
        pileup_file = open(pileup_filename, "r")
        for line in pileup_file:
            line = line.strip()
            elements = line.split("\t")
            tuple = (elements[1], elements[2])
            list.append(tuple)


        match = re.match(r"^(.*)_(\d+).*", pileup_filename)
        minifile_name = match.group(1)
        new_file = open(minifile_name + "_normalized_coverage.bed", "w")
        # Store coord range from "sliding_window_annotation file"
        window_annotation_file = open("windows_51_20_strandp.bed", "r")
        for line in window_annotation_file:
            line = line.strip()
            if line.startswith("initial_coord"):
                continue
            elements = line.split("\t")
            initial_coord = int(elements[0])
            final_coord = int(elements[1])
            window_range = range(initial_coord, final_coord+1)
            sum_coverage = 0
            iteration = 0
            for tuple in list:
                if int(tuple[0]) in window_range:
                    iteration += 1
                    perbase_coverage = float(tuple[1])
                    sum_coverage = sum_coverage + perbase_coverage
                if iteration == 51:
                    window_size = (final_coord-initial_coord)+1
                    mean_coverage = sum_coverage / window_size
                    normalized_mean_coverage = mean_coverage / sum_ribosomal_cov
                    del list[0:18]
                    break
            print(line + "\t" + str(normalized_mean_coverage))
            new_file.write(line + "\t" + str(normalized_mean_coverage) + "\n")

calculate_normalization()


def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]


def merging_normalized_coverages():
    file_names = glob.glob("*normalized_coverage.bed")
    merged_file = open("merged_coverage.bed", "w")
    coverage_dict = dict()
    minifile_name_list = []
    common_information_list = []
    for file_name in file_names:
        file = open(file_name, "r")
        match = re.match(r"^(.*)_normalized", file_name)
        minifile_name = match.group(1)
        minifile_name_list.append(minifile_name)
        coverage_list = []
        for line in file:
            line = line.strip()
            elements = line.split("\t")
            normalized_coverage = elements[10]
            common_information = elements[0] + "\t" + elements[1] + "\t" +elements[2] + "\t" +elements[3] + "\t" +elements[4] + "\t" +elements[5] + "\t" +elements[6] + "\t" +elements[7] + "\t" + elements[8] + "\t" + elements[9]
            coverage_list.append(normalized_coverage)
            if file_name is file_names[-1]:
                common_information_list.append(common_information)
        coverage_dict[minifile_name] = coverage_list

    minifile_name_list.sort(key=natural_keys)

    # Writing data
    merged_file.write("\t\t\t\t\t\t\t\t\t\t")
    for minifile_name in minifile_name_list:
        merged_file.write(minifile_name)
        if minifile_name is not minifile_name_list[-1]:
            merged_file.write("\t")
    merged_file.write("\n")
    for common_information in common_information_list:
        merged_file.write(common_information + "\t")
        for minifile_name in minifile_name_list:
            merged_file.write(coverage_dict[minifile_name][0])
            if minifile_name is not minifile_name_list[-1]:
                merged_file.write("\t")
            coverage_dict[minifile_name].remove(coverage_dict[minifile_name][0])
        merged_file.write("\n")

merging_normalized_coverages()


def fitting_regression():
    file = open("merged_coverage.bed", "r")
    new_file = open("half-lives.bed", "w")
    new_file.write("initial_coord\tfinal_coord\tstrand\tgenes\tattributes\tnumber_TSS\toperons\tTSS_position\tCDS_initial_position\tCDS_final_position\tk(decay)\thalf_live\tr_squared\n")

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
            common_information = elements[0] + "\t" + elements[1] + "\t" +elements[2] + "\t" +elements[3] + "\t" +elements[4] + "\t" +elements[5] + "\t" +elements[6] + "\t" +elements[7] + "\t" + elements[8] + "\t" + elements[9]
            y_list = elements[10:]
            y_list = list(map(float, y_list))
            for number in y_list:
                if number == 0:
                    index = y_list.index(number)
                    y_list[index] = 0.000000000000000000000000002

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

            # Calculate half-lives:
            half_live = np.log(2) / (-k_decay)
            #print(common_information + "\nK decay=" + str(-k_decay) + "\tc=" + str(c) + "\tHalf-live=" + str(half_live) + "r-squared=" + str(R2))

            # Visalizing data:
            y = np.exp(k_decay*x_data) * np.exp(c)
            plt.plot(x_data, y_data, "o")
            plt.plot(x_data, y)
            plt.title(elements[3] + "   (" + elements[0] + "-" + elements[1] + ")\n" + "operons: " + elements[6] + "  " + "number TTS: " + elements[5] + "\n" +  "half-lfe: " + str(half_live) + "   r-squared:" + str(R2))
            #plt.show()

            new_file.write(common_information + "\t" + str(-k_decay) + "\t" + str(half_live) + "\t" + str(R2) + "\n")
fitting_regression()
