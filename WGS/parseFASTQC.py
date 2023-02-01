

# ---
# --- About: 
# --- 
# ---   Programm that parses all fastqc.html files in curren directory and create a csv file with results
# ---   Rows are samples 
# ---   Collumns are Summary statistics - Basic Statistics, Per base sequence quality, per sequence quality scores, per base sequence content, per sequence GC content, per base N content, sequence length distribution, sequence duplication levels, overrepresented sequences, adapter content
# ---   Collumns are also Total Sequences, Sequences flagged as poor quality, %GC
# --- IN: 
# ---   path to fastqc.html files
# ---   result csv file name/path
# --- OUT:
# ---   csv file with results at CSV file path
# --- 
# --- Author: Edgars Liepa edgars.liepa@biomed.lu.lv
# --- Date: 18.1.23
# --- Modified: 

import os
import re
import csv
import pandas as pd
import sys
import argparse

def parse_fastqc(args, parser):

    # Create a list of all files in current directory
    if not args.path:
        files = os.listdir(os.getcwd())
    else:
        files = os.listdir(args.path)

    # Create a list of all fastqc.html files

    fastqc = []
    for file in files:
        if file.endswith(".html"):
            fastqc.append(file)

    # Create a list of all samples remowing _fastqc.html
    samples = []
    for file in fastqc:
        sample = re.sub("_fastqc.html", "", file)
        samples.append(sample)

    # Create column names for csv file with Summary statistics
    statistics = ["Sample Name", "Total Sequences", "Sequences flagged as poor quality", "%GC", "Basic Statistics", "Per base sequence quality", "Per sequence quality scores", "Per base sequence content", "Per sequence GC content", "Per base N content", "Sequence Length Distribution", "Sequence Duplication Levels", "Overrepresented sequences", "Adapter Content"]

    from bs4 import BeautifulSoup

    # Parse .html files and extract for every statistics if it was failed, passed or warning
    # create a row with sample name and results for each sample
    # find a HTML tag  alt="[status]" and extract status
    # add status to collumn with statistics

    df = []
    # Create emtpt data frame with column names
    df = pd.DataFrame(columns=statistics)

    for file in fastqc:
        print("Process: ", file)
        with open(file, "r") as f:
            # add sample name to df Sample Name collumn 
            sample = re.sub("_fastqc.html", "", file)
            # create new row with sample name in df and use concat to add it to df
            df = pd.concat([df, pd.DataFrame([sample], columns=["Sample Name"])], ignore_index=True)

            # parse HTML using html parser
            soup = BeautifulSoup(f, "html.parser")
            # find in class summary find all href tags
            
            for link in soup.find(class_="summary").find_all("li"):
                
                rez = link.find("img").get("alt")
                # print(rez)
                
                # find text in tag a and add it to tag
                tag = link.find("a").text            
                # print(tag)
                
                # add rez to tag column in df   
                df.loc[df["Sample Name"] == sample,tag] = rez
                
                
            # find table within soup and extract Measure and Value columns
            for row in soup.find("table").find_all("tr"):
                
                # find all td tags
                value = row.find_all("td")

                # if value non empty
                if value:
                    # print(value)
                    # print(value[1].text)
                    # add value to df row with Sample name column value sample and row value[0]
                    df.loc[df["Sample Name"] == sample, value[0].text] = value[1].text


    # Save dataframe df to csv file
    if not args.csv:
        df.to_csv("resultsFASTQC.csv", index=False)
    else:
        df.to_csv(args.csv + "/resultsFASTQC.csv", index=False)

# create main function
def main():
    # get two arguments from command line
    # 1. path to fastqc.html file folder
    # 2. path to csv file

    parser = argparse.ArgumentParser(description = 'Programm that parses all fastqc.html files in a directory and create a csv file with results',
                                    epilog = 'If no arguments are given, current directory is used and csv file is created in current directory with name resultsFASTQC.csv')
    parser.add_argument("--path", help="path to fastqc.html file directory")
    parser.add_argument("--csv", help="path to csv file")
    args = parser.parse_args()

    # # if arguments are empty print help
    # if not args.path:
    #     parser.print_help()
    #     sys.exit(1) 
    
    # if not args.csv:
    #     parser.print_help()
    #     sys.exit(1)
    
    print(args.path)
    print(args.csv)
    parse_fastqc(args, parser)

    if not args.csv:
        print("Done \nResults are saved in {} resultsFASTQC.csv" .format(os.getcwd()))
    else:
        print("Done \nResults are saved in {}/resultsFASTQC.csv file" .format(args.path))


if __name__ == "__main__":
    main()