# Author: Louis-Jacques Ruel (louis.jacques.ruel@gmail.com)
# Last modification: 2023-11-06

# Executing SigProfilerExtractor to get SBS signatures of all patients.
# Reference : https://github.com/AlexandrovLab/SigProfilerExtractor
from SigProfilerExtractor import sigpro as sig

def main_function():
    # Main function with default values, except cpu usage and maximum_signatures.
    sig.sigProfilerExtractor("vcf", "NSLC_SBS", "92_patients_NSLC_filtered_VCFS", "GRCh37", cpu=4, minimum_signatures=1, maximum_signatures=8)

if __name__=="__main__":
   main_function()


