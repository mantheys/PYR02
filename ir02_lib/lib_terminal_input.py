import numpy as np

def terminal_input():
    run_list = []; ch_list = []
    run_list_string = input("Input run list (use space between numbers): ")
    ch_list_string = input("Input channel list (use space between numbers): ")
    run_list_split = run_list_string.split(" ")
    ch_list_split = ch_list_string.split(" ")

    for i in range(np.size(run_list_split)):
        run_list.append(int(run_list_split[i])) 

    for i in range(np.size(ch_list_split)):
        ch_list.append(int(ch_list_split[i]))
    
    return run_list, ch_list