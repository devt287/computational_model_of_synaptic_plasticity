
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import os
import pandas as pd
from IPython.display import display


#import the txt file to the dictionary
def import_data(filepath):
    # Initialize an empty dictionary
    results = {}

    # Open the file and read line by line
    with open(filepath, 'r') as file:
        for line_number, line in enumerate(file, 1):  # Added line_number for debugging
            # Split the line to isolate the part with the numbers
            try:
                # Assuming each line follows the format: "G_VDCC=X.XX, G_NMDA=Y.YY, [list of numbers]"
                g_vdcc_part, g_nmda_part, numbers_part = line.split(',', 2)
                # Extract G_VDCC value
                g_vdcc = float(g_vdcc_part.split('=')[1].strip())
                # Extract numbers from the last part
                numbers_str = numbers_part.strip()[1:-1]  # Strip the enclosing brackets and whitespace
                data_list = [int(num.strip()) for num in numbers_str.split(',') if num.strip()]
            except ValueError as e:
                print(f"Error parsing on line {line_number}: {e}")
                continue

            # Add to the dictionary
            results[g_vdcc] = data_list

    return results

#sort the txt
def export_data(data_dict, filepath):
    """
    Exports the data from a dictionary to a file.

    Args:
    data_dict (dict): The dictionary containing the data to export.
    filepath (str): The path to the file where the data should be written.
    """
    sorted_dict = sort_dict_by_keys(data_dict)

    with open(filepath, 'w') as file:
        for g_vdcc, data_list in sorted_dict.items():
            g_vdcc_formatted = f"{g_vdcc:.2f}"
            g_nmda = round(g_vdcc * 3, 2)
            g_nmda_formatted = f"{g_nmda:.2f}"
            numbers_str = ', '.join(map(str, data_list))
            line = f"G_VDCC={g_vdcc_formatted}, G_NMDA={g_nmda_formatted}, [{numbers_str}]\n"
            file.write(line)

def sort_dict_by_keys(data_dict):
    """
    Sorts a dictionary by its keys and returns a new dictionary.

    Args:
    data_dict (dict): The dictionary to sort.

    Returns:
    dict: A new dictionary sorted by keys.
    """
    sorted_dict = {k: data_dict[k] for k in sorted(data_dict)}
    return sorted_dict


def zero_wing(data_dict):
    # List to hold keys that meet the criteria
    target_keys = []
    
    # Iterate over each key and list in the dictionary
    for key, values in data_dict.items():
        # Check if the list starts and ends with zero and is not empty
        if values and values[0] == 0 and values[-1] == 0:
            target_keys.append(key)
    
    return target_keys

def DP(data_dict):
    # List to hold keys that meet all criteria
    target_keys = []

    # Iterate over each key and list in the dictionary
    for key, values in data_dict.items():
        # Check if the list starts and ends with zero
        if values and values[0] == 0 and values[-1] == 0:
            # Check the sequence of -1 before 1 after the initial zero
            seen_minus_one = False
            seen_one = False
            valid_sequence = True  # Assume the sequence is valid initially
            
            # Iterate through elements in the list except the first and last
            for value in values[1:-1]:
                if value == -1:
                    if seen_one:  # If 1 has been seen and -1 comes again
                        valid_sequence = False
                        break
                    seen_minus_one = True
                elif value == 1:
                    if seen_minus_one:  # Once we've seen -1, 1 can appear
                        seen_one = True
                    else:  # If 1 comes before any -1
                        valid_sequence = False
                        break
            
            # Only add key if the sequence is valid
            if valid_sequence and seen_minus_one and seen_one:
                target_keys.append(key)

    return target_keys


def DPD(data_dict):
    # List to hold keys that meet all criteria
    target_keys = []

    # Iterate over each key and list in the dictionary
    for key, values in data_dict.items():
        # Check if the list starts and ends with zero
        if values and values[0] == 0 and values[-1] == 0:
            # Initiate state tracking, where we haven't seen -1 or 1 yet
            seen_minus_one_before_one = False
            seen_one = False
            seen_minus_one_after_one = False
            
            # Iterate through elements in the list except the first and last
            for value in values[1:-1]:
                if value == -1:
                    if not seen_one:
                        seen_minus_one_before_one = True  # First occurrence of -1 before any 1
                    else:
                        seen_minus_one_after_one = True  # Occurrence of -1 after 1
                elif value == 1:
                    if seen_minus_one_before_one:  # We have already seen -1
                        seen_one = True
                    else:
                        # If 1 appears before any -1, break and do not add to list
                        break
            
            # Check final condition to add the key
            if seen_minus_one_before_one and seen_one and seen_minus_one_after_one:
                target_keys.append(key)

    return target_keys

def DPD_prime(data_dict):
    # List to hold keys that meet all criteria
    target_keys = []

    # Iterate over each key and list in the dictionary
    for key, values in data_dict.items():
        # Check if the list starts and ends with -1
        if values and values[0] == -1 and values[-1] == -1:
            # Find the first index of 1 after the initial -1s
            try:
                first_one_index = values.index(1)
                # Find the last index of 1 before the final -1s
                last_one_index = len(values) - 1 - values[::-1].index(1)
            except ValueError:
                # If there is no 1 in the list, skip this list
                continue

            # Check that all elements before the first 1 are -1
            if all(value == -1 for value in values[:first_one_index]):
                # Check that all elements after the last 1 are -1
                if all(value == -1 for value in values[last_one_index + 1:]):
                    # Check that all elements between the first and last 1 are 1
                    if all(value == 1 for value in values[first_one_index:last_one_index + 1]):
                        target_keys.append(key)

    return target_keys


def all_zero(data_dict):
    # List to hold keys where all entries in the list are zero
    zero_keys = []

    # Iterate over each key and its associated list in the dictionary
    for key, values in data_dict.items():
        # Check if all items in the list are zero using the all() function
        if all(value == 0 for value in values):
            zero_keys.append(key)

    return zero_keys

def PD(data_dict):
    # List to hold keys that meet all criteria
    target_keys = []

    # Iterate over each key and list in the dictionary
    for key, values in data_dict.items():
        # Check if the list starts and ends with zero
        if values and values[0] == 0 and values[-1] == 0:
            # Check the sequence of 1 before -1 after the initial zero
            seen_one = False
            seen_minus_one = False
            valid_sequence = True  # Assume the sequence is valid initially
            
            # Iterate through elements in the list except the first and last
            for value in values[1:-1]:
                if value == 1:
                    if seen_minus_one:  # If -1 has been seen and 1 comes again
                        valid_sequence = False
                        break
                    seen_one = True
                elif value == -1:
                    if seen_one:  # Once we've seen 1, -1 can appear
                        seen_minus_one = True
                    else:  # If -1 comes before any 1
                        valid_sequence = False
                        break
            
            # Only add key if the sequence is valid
            if valid_sequence and seen_one and seen_minus_one:
                target_keys.append(key)

    return target_keys




def all_one(data_dict):
    # List to hold keys where all entries in the list are one
    one_keys = []

    # Iterate over each key and its associated list in the dictionary
    for key, values in data_dict.items():
        # Check if all items in the list are one using the all() function
        if all(value == 1 for value in values):
            one_keys.append(key)

    return one_keys


def all_negative_one(data_dict):
    # List to hold keys where all entries in the list are negative one
    negative_one_keys = []

    # Iterate over each key and its associated list in the dictionary
    for key, values in data_dict.items():
        # Check if all items in the list are negative one using the all() function
        if all(value == -1 for value in values):
            negative_one_keys.append(key)

    return negative_one_keys


def P(data_dict):
    # List to hold keys where the list starts and ends with zero and has at least one '1' in the middle
    target_keys = []

    # Iterate over each key and its associated list in the dictionary
    for key, values in data_dict.items():
        # Check if the list starts and ends with zero
        if values and values[0] == 0 and values[-1] == 0:
            # Check if there is at least one '1' in the middle
            if 1 in values[1:-1]:  # Slicing to exclude the first and last element
                target_keys.append(key)

    return target_keys


def D(data_dict):
    # List to hold keys where the list starts and ends with zero and has any sequence of continuous '-1's flanked by '0's in the middle
    target_keys = []

    # Iterate over each key and its associated list in the dictionary
    for key, values in data_dict.items():
        # Check if the list starts and ends with zero
        if values and values[0] == 0 and values[-1] == 0:
            # Check for sequences of -1s and ensure no 1s are present
            valid_sequence_found = False
            i = 1  # Start from the first middle element

            while i < len(values) - 1:
                if values[i] == -1:
                    # Start of a potential -1 sequence
                    start = i
                    while i < len(values) - 1 and values[i] == -1:
                        i += 1
                    # Ensure the sequence of -1 is flanked by 0s
                    if (start == 1 or values[start - 1] == 0) and (i == len(values) - 1 or values[i] == 0):
                        valid_sequence_found = True
                elif values[i] == 1:
                    # If a 1 is found, invalidate the sequence
                    valid_sequence_found = False
                    break
                else:
                    i += 1

            if valid_sequence_found:
                target_keys.append(key)

    return target_keys


def P(data_dict):
    # List to hold keys where the list starts and ends with zero and has any sequence of continuous '1's flanked by '0's in the middle
    target_keys = []

    # Iterate over each key and its associated list in the dictionary
    for key, values in data_dict.items():
        # Check if the list starts and ends with zero
        if values and values[0] == 0 and values[-1] == 0:
            # Variable to track if a valid sequence of 1s has been found
            valid_sequence_found = False
            i = 1  # Start from the first middle element

            while i < len(values) - 1:
                if values[i] == 1:
                    # Start of a potential '1' sequence
                    start = i
                    while i < len(values) - 1 and values[i] == 1:
                        i += 1
                    # Ensure the sequence of '1' is flanked by '0s'
                    if (start == 1 or values[start - 1] == 0) and (i == len(values) - 1 or values[i] == 0):
                        valid_sequence_found = True
                        break
                else:
                    i += 1

            if valid_sequence_found:
                target_keys.append(key)

    return target_keys


def have_overlap(*lists):
    # Convert the first list to a set
    common_elements = set(lists[0])
    
    # Iterate over the remaining lists and update the set with intersection
    for lst in lists[1:]:
        common_elements.intersection_update(lst)
        
        # If at any point the set of common elements becomes empty, return False
        if not common_elements:
            return False
    
    # If common elements exist after checking all lists, return True
    return bool(common_elements)

def get_info(data_dict):
    DP_keys = DP(data_dict)
    PD_keys = PD(data_dict)
    DPD_keys = DPD(data_dict)
    DPD_prime_keys = DPD_prime(data_dict)
    P_keys = P(data_dict)
    D_keys = D(data_dict)
    all_zero_keys = all_zero(data_dict)
    all_one_keys = all_one(data_dict)
    all_negative_one_keys = all_negative_one(data_dict)
    return [D_keys,P_keys,DPD_keys,DPD_prime_keys,all_zero_keys,all_one_keys,all_negative_one_keys,DP_keys,PD_keys]



def visualize_lists(lists, list_names=["D", "P", "DPD", "DPD'", "all_basal", "P'", "D'","DP","PD"], colors=None):
    """
    Visualizes multiple lists on a single plot with each list represented by vertical bands of a unique color,
    including an explicit 'Other' category in the legend for non-covered areas, with enhanced x-axis grid.

    Args:
    lists (list of lists): Each sublist contains numerical values to be plotted.
    list_names (list of str): Names for each list for the legend, automatically adjusted to remove '_list'.
    colors (list of str, optional): Colors for each list. If None, colors will be chosen automatically.
    """
    # Check if colors are provided, otherwise generate them
    if colors is None:
        # Generate a colormap
        cmap = plt.get_cmap('tab20')
        colors = [cmap(i) for i in np.linspace(0, 1, len(lists))]

    # Initialize the plot
    fig, ax = plt.subplots()

    # Plot each list with its respective color
    for data, color, name in zip(lists, colors, list_names):
        for value in data:
            ax.axvline(x=value, color=color, linewidth=2, label=name if name not in ax.get_legend_handles_labels()[1] else "")

    # Adding labels and legend
    ax.set_xlim(0, 10)  # Set the range from 0 to 10
    ax.set_xlabel('G_VDCC')
    ax.set_ylabel('Presence')
    ax.set_title('Final State')

    # Set x-axis ticks and labels
    ax.set_xticks(np.arange(0, 10.1, 0.5))
    ax.set_xticklabels([f"{val:.1f}" for val in np.arange(0, 10.1, 0.5)])

    # Enable grid for better precision
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.minorticks_on()

    # Customizing the legend to omit 'list' and capitalize
    handles, labels = plt.gca().get_legend_handles_labels()

    # Create a custom legend item for 'Other'
    other_patch = mpatches.Patch(color='white', label='Other')
    handles.append(other_patch)  # Append the custom patch
    labels.append('Other')  # Append the label for 'Other'

    # Remove duplicates and create final legend
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), [label.replace('_list', '').upper() for label in by_label.keys()])

    
    # Show the plot
    plt.show()


def visualize_list_lengths(lists, list_names, colors=None):
    """
    Generates a bar chart where each bar represents the length of a list.

    Args:
    lists (list of lists): Each sublist whose length is to be plotted.
    list_names (list of str): Names for each list to be used in the legend.
    colors (list of str, optioinal): Colors for each bar. If None, colors will be chosen automatically.
    """
    # Extract the lengths of each list
    list_lengths = [len(lst) for lst in lists]
    
    # Check if colors are provided, otherwise generate them
    if colors is None:
        # Generate a colormap
        cmap = plt.get_cmap('tab20')
        colors = [cmap(i) for i in np.linspace(0, 1, len(lists))]

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Create a bar chart
    bars = ax.bar(list_names, list_lengths, color=colors)

    # Set the title and labels

    ax.set_title('Bar Chart of Each State')

    # Set y-axis to show exact numbers
    ax.set_yticks(list_lengths)
    ax.set_yticklabels([str(length) for length in list_lengths])

    # Customizing the legend to simplify names
    handles, labels = bars, [name.replace('_list', '').upper() for name in list_names]
    ax.legend(handles, labels)

    # Display the plot
    plt.show()


def display_grid_form(data_dict, keys, file_name=None):
    header = [-500, -400, -300, -200, -150, -100, -75, -50, -25, -10, 0, 10, 25, 50, 75, 100, 200, 300, 400, 500]
    data = []

    for key in keys:
        if key in data_dict:
            data.append([key] + data_dict[key])
        else:
            print(f"Key {key} not found in the dictionary.")
            data.append([key] + [None] * len(header))

    # Create DataFrame with header
    df = pd.DataFrame(data, columns=["G_VDCC"] + header)
    
    # Set display options to show all columns
    pd.set_option('display.max_columns', None)
    
    # Display the DataFrame
    display(df)
    
    # Save DataFrame to a CSV file if file_name is provided
    if file_name:
        df.to_csv(file_name, index=False)
        print(f"Data saved to {file_name}")



def visualize_multiple_list_groups(y_labels, list_groups, names, list_names=["D", "P", "DPD", "DPD'", "all_basal", "P'", "D'", "DP", "PD"], colors=None):
    """
    Visualizes multiple groups of lists on a single plot with each group represented by vertical bands of unique colors,
    with enhanced x-axis grid.

    Args:
    y_labels (list of str): Labels for each group to be displayed on the y-axis.
    list_groups (list of lists of lists): A list containing lists of lists of numerical values to be plotted.
    list_names (list of str): Names for each list for the legend, automatically adjusted to remove '_list'.
    colors (list of str, optional): Colors for each list. If None, colors will be chosen automatically.
    """
    num_groups = len(list_groups)
    if num_groups == 0:
        raise ValueError("At least one list group must be provided.")
    if len(y_labels) != num_groups:
        raise ValueError("The number of y_labels must match the number of list groups.")

    # Check if colors are provided, otherwise generate them
    if colors is None:
        # Generate a colormap
        cmap = plt.get_cmap('tab20')
        colors = [cmap(i) for i in np.linspace(0, 1, len(list_names))]

    # Initialize the plot
    fig, ax = plt.subplots()

    # Calculate the height of each group section on the y-axis
    section_height = 1 / num_groups

    # Find the minimum of the maximum values of G_VDCC for each row
    max_g_vdcc_per_group = [max(max(lists, default=[0]), default=0) for lists in list_groups]
    min_of_max_g_vdcc = min(max_g_vdcc_per_group)

    # Plot each group with its respective color
    for group_index, lists in enumerate(list_groups):
        for data, color, name in zip(lists, colors, list_names):
            for value in data:
                y_min = group_index * section_height
                y_max = (group_index + 1) * section_height
                ax.axvline(x=value, ymin=y_min, ymax=y_max, color=color, linewidth=2, label=name if name not in ax.get_legend_handles_labels()[1] else "")

    # Adding labels and legend
    ax.set_xlim(0, min_of_max_g_vdcc)  # Set the range from 0 to the minimum of max values of G_VDCC
    ax.set_xlabel('G_VDCC')
    ax.set_ylabel(names)
    ax.set_title('Final State of different ' + names)

    # Set x-axis ticks and labels
    ax.set_xticks(np.arange(0, min_of_max_g_vdcc + 0.5, 0.5))
    ax.set_xticklabels([f"{val:.1f}" for val in np.arange(0, min_of_max_g_vdcc + 0.5, 0.5)])

    # Set y-axis labels
    y_ticks = [(i + 0.5) * section_height for i in range(num_groups)]
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)

    # Enable grid for better precision
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.minorticks_on()

    # Adding horizontal lines to separate groups
    for i in range(num_groups + 1):
        ax.axhline(y=i * section_height, color='black', linestyle='--', linewidth=0.8)

    # Customizing the legend to omit 'list' and capitalize
    handles, labels = plt.gca().get_legend_handles_labels()

    # Create a custom legend item for 'Other'
    other_patch = mpatches.Patch(color='white', label='Other')
    handles.append(other_patch)  # Append the custom patch
    labels.append('Other')  # Append the label for 'Other'

    # Remove duplicates and create final legend
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), [label.replace('_list', '').upper() for label in by_label.keys()], loc='upper right')

    # Show the plot
    plt.show()


def visualize_each_element(info_list, element_list, name_list):
    """
    Visualizes counts of each element in the provided info lists on a line plot.

    Args:
    info_list (list of lists): Each sublist contains multiple lists representing different elements.
    element_list (list of str): Names for each element to be displayed on the x-axis.
    name_list (list of str): Names to be used in the legend to denote each list in info_list.
    """
    num_elements = len(element_list)
    num_info = len(info_list)
    
    # Create a colormap
    cmap = plt.get_cmap('tab20')
    colors = [cmap(i) for i in np.linspace(0, 1, num_info)]
    
    # Initialize the plot with larger size
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Iterate over each info list and plot it
    for idx, info in enumerate(info_list):
        counts = [len(info[i]) if i < len(info) else 0 for i in range(num_elements)]
        ax.plot(element_list, counts, label=name_list[idx], color=colors[idx], marker='o')
    
    # Set the title and labels
    ax.set_title(f'Counts of Each Element in {num_info} Lists')
    ax.set_xlabel('Elements')
    ax.set_ylabel('Count')
    
    # Adding legend
    ax.legend(loc='upper right')
    
    # Enable grid for better precision
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.minorticks_on()
    
    # Display the plot
    plt.show()


def display_DP(data_dict_list, name_list, header, column_names, file_name=None):
    combined_data = []

    for data_dict, name in zip(data_dict_list, name_list):
        dp_keys = DP(data_dict)
        if len(dp_keys) > 0:
            for key in dp_keys:
                combined_data.append([name, key] + data_dict[key])

    # Create DataFrame with provided column names
    df = pd.DataFrame(combined_data, columns=column_names + header)
    
    # Set display options to show all columns
    pd.set_option('display.max_columns', None)
    
    # Display the DataFrame
    display(df)
    
    # Save DataFrame to an Excel file if file_name is provided
    if file_name:
        df.to_excel(file_name, index=False, engine='xlwt')
        print(f"Data saved to {file_name}")

def display_DPD(data_dict_list, name_list, header, column_names, file_name=None):
    combined_data = []

    for data_dict, name in zip(data_dict_list, name_list):
        dp_keys = DPD(data_dict)
        if len(dp_keys) > 0:
            for key in dp_keys:
                combined_data.append([name, key] + data_dict[key])

    # Create DataFrame with provided column names
    df = pd.DataFrame(combined_data, columns=column_names + header)
    
    # Set display options to show all columns
    pd.set_option('display.max_columns', None)
    
    # Display the DataFrame
    display(df)
    
    # Save DataFrame to an Excel file if file_name is provided
    if file_name:
        df.to_excel(file_name, index=False, engine='xlwt')
        print(f"Data saved to {file_name}")


def visualize_multi_list_lengths(list_of_lists, list_names, name_list, save_address=os.getcwd()):
    """
    Generates line charts for each sublist in the list of lists, where each line represents the length of the lists.

    Args:
    list_of_lists (list of lists of lists): Each entry in the outer list is a list of lists to be plotted.
    list_names (list of str): Names for each element in the x-axis.
    name_list (list of str): Names for each group to be used in the legend.
    save_address (str): The directory where the plots will be saved. Defaults to the current directory.
    """
    # Check if colors are provided, otherwise generate them
    cmap = plt.get_cmap('tab20')
    colors = [cmap(i) for i in np.linspace(0, 1, len(list_of_lists))]

    # Iterate through each list and create a separate plot
    for idx, (lists, name) in enumerate(zip(list_of_lists, name_list)):
        list_lengths = [len(lst) for lst in lists]
        
        # Create a figure and axis
        fig, ax = plt.subplots(figsize=(14, 10))  # Increased size for better visibility

        # Plot the lengths
        ax.plot(list_names, list_lengths, label=name, color=colors[idx % len(colors)], marker='o')

        # Set the title and labels
        ax.set_title(f'Line Chart of Each State for {name}')
        ax.set_xlabel('States')
        ax.set_ylabel('Count')
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)

        # Add a legend at the top right corner
        ax.legend(loc='upper right')

        # Save the plot to the specified directory
        plt.savefig(os.path.join(save_address, f'name_list[idx].png'))
        plt.close()


def DPD_pattern(a, a_info, save_path, show="true"):
    """
    Generate a horizontal heatmap with adjusted entry dimensions based on the DPD cases from a_info.
    The x-axis corresponds to delta_t_values, and the y-axis corresponds to sorted DPD cases' G_VDCC values.

    Parameters:
    - a: Dictionary from `import_data`, containing G_VDCC values as keys and lists as values.
    - a_info: List from the `get_info` function, containing information about DPD cases and other parameters.
    - save_path: Path to save the generated heatmap.
    - show: If "true", display the plot. If not "true", skip displaying it.
    """
    # Define delta_t_values
    delta_t_values = [-500, -400, -300, -200, -150, -100, -75, -50, -25, -10, 0, 10, 25, 50, 75, 100, 200, 300, 400, 500]

    # Extract DPD keys from a_info
    DPD_keys = a_info[2]  # Assuming a_info[2] contains DPD keys

    # Sort G_VDCC values (DPD keys) for x-axis
    G_VDCC_values = sorted([dpd for dpd in DPD_keys if dpd in a])

    # Initialize the matrix (x-axis: G_VDCC cases, y-axis: delta_t_values)
    matrix = np.zeros((len(delta_t_values), len(G_VDCC_values)), dtype=int)

    # Populate the matrix
    for x_idx, g_vdcc in enumerate(G_VDCC_values):
        if g_vdcc in a:
            column_values = a[g_vdcc]
            matrix[:, x_idx] = column_values  # Assign the column values from a[g_vdcc]

    # Create the heatmap
    fig, ax = plt.subplots(figsize=(16, 8))  # Adjust the overall figure size
    cax = ax.matshow(matrix.T, cmap='RdBu', vmin=-1, vmax=1, aspect=0.5)  # Aspect=2 for 1:2 height:width ratio

    # Add colorbar with labels
    cbar = plt.colorbar(cax, ticks=[-1, 0, 1])
    cbar.ax.set_yticklabels(['LTD', 'Basal', 'LTP'])

    # Set axis labels and ticks
    ax.set_xticks(range(len(delta_t_values)))
    ax.set_xticklabels(delta_t_values, rotation=0)  # Keep delta_t values horizontal
    ax.set_yticks(range(len(G_VDCC_values)))
    ax.set_yticklabels([f"{g_vdcc:.2f}" for g_vdcc in G_VDCC_values])
    ax.set_xlabel('Delta t values')
    ax.set_ylabel('G_VDCC')

    # Add title
    plt.title('DPD Heatmap plot', pad=20)
    plt.tight_layout()

    # Show the plot only if show == "true"
    if show == "true":
        plt.show()

    # Save the plot
    plt.savefig(save_path)
    plt.close()
    print(f"Saved heatmap to {save_path}")


def DPD_combo(a_list, a_info_list, path):
    """
    Generate heatmaps for each DPD pattern in the a_list and save them to the given path.

    Parameters:
    - a_list: List of data dictionaries.
    - a_info_list: List of info data corresponding to a_list.
    - path: Directory path to save the heatmaps.
    """
    # Ensure the directory exists
    if not os.path.exists(path):
        os.makedirs(path)

    # Loop through the lists and generate heatmaps
    for i, (data, info) in enumerate(zip(a_list, a_info_list)):
        # Use the provided path to save the heatmap
        save_path = os.path.join(path, f"heatmap_{i + 1}.png")
        DPD_pattern(data, info, save_path,"false")



