from ast import literal_eval # for converting strings to integers
import re
from typing import List
import csv
import os



# run class summarizing each output block
class Run:
    def __init__(self, call, solved, time):
        # e.g. stdin = "ffpo-transport-p17"
        stdin = re.split("-", call)
        self.planner = stdin[0]                                             # "ffpo"
        self.domain = ""
        if len(stdin)==3:
            self.domain = stdin[1]
        else:
            self.domain = self.domain+stdin[1]
            for i in range(2,len(stdin)-1):
                self.domain = self.domain+"-"+stdin[i]
        self.instance = -2
        to_eval = stdin[len(stdin)-1][slice(1,len(stdin[2]))]
        if to_eval[0] == '0':
            self.instance = literal_eval( to_eval[slice(1,len(to_eval))] )      # "17"
        else:
            self.instance = literal_eval(to_eval)
        self.solved = solved
        self.time = time


runs = []
path = "./runs_2.txt"
with open(path, "r") as file:
    content = file.read()
    lines = content.split("\n")
    numlines = len(lines)
    info_pad_len = len("INFO     search stdin: ")
    running_str = "INFO     Running s"
    search_exit_str = "search e"

    i=0
    count = 0
    while i < numlines:
        line = lines[i]
        if line[slice(0,18)] == running_str:                # if at start of run
            i = i+1
            call = lines[i][slice(info_pad_len,-4)]         # get the call used e.g. "ffpo-transport-p17"

            while i < numlines:                             # now look for exit code
                i = i+1
                line = lines[i]
                if (line[slice(17)] == "search exit code:"): # bad run found
                    if (line[-1] != '0'):
                        i = i+1
                        break
                    elif line == "search exit code: 0":     # good run found, look for search time
                        t_line = lines[i-4]
                        colon_ind = t_line.find(':')
                        runtime = literal_eval(t_line[slice(colon_ind+2, -1)])
                        run = Run(call, True, runtime)
                        runs.append(run)
                        count = count+1
                        i = i+1
                        break
        else:
            i = i+1
            


def unique_planners(runs: List[Run]):
    unique_list = []
    for run in runs:
        if run.planner not in unique_list:
            unique_list.append(run.planner)
    return unique_list


def unique_domains(runs: List[Run]):
    unique_list = []
    for run in runs:
        if run.domain not in unique_list:
            unique_list.append(run.domain)
    return unique_list


def unique_instances(runs: List[Run]):
    unique_list = []
    for run in runs:
        if run.instance not in unique_list:
            unique_list.append(run.instance)
    return unique_list

planners = unique_planners(runs)
domains = unique_domains(runs)
instances = unique_instances(runs)

# for each (domain, planner) pair, create a csv file and write the points to it
for planner in planners:
    for dom in domains:
        dom_plan_runs = []
        for run in runs:
            if (run.planner == planner) and (run.domain == dom):
                dom_plan_runs.append(run)

        file_name = "/planner_data/"+planner+"-"+dom+".csv"
        with open(os.getcwd()+file_name, "w+", newline='') as file:
            writer = csv.writer(file, delimiter=",")
            writer.writerow(["instance", "runtime"])
            for run in dom_plan_runs:
                row = [run.instance, run.time]
                writer.writerow(row)

        file.close()