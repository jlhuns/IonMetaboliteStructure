import file_paths as FILE_PATHS
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

def load_and_categorize_data(file_path):
    # Load the CSV data
    data = pd.read_csv(file_path)
    
    # Separate into Ion and Metabolite groups based on the Ligand_Type column
    ion_data = data[data['Ligand_Type'] == 'Ion'].copy()
    metabolite_data = data[data['Ligand_Type'] == 'Metabolite'].copy()
    
    return ion_data, metabolite_data

def compute_conservation_scores(data):
    # Map each symbol to its corresponding numerical value.
    symbol_to_value = {
        '*': 1.0,  # Highly conserved
        ':': 0.8,  # Moderately conserved
        '.': 0.5,  # Low conservation
        '-': 0.0   # No conservation (gap)
    }
    
    # Convert the conservation score strings to lists of numerical values.
    data['conservationScore'] = data['conservationScore'].apply(
        lambda score: [symbol_to_value[char] for char in score if char in symbol_to_value]
    )
    return data

def compute_average_per_entry(data):
    # Compute the average score for each entry.
    data['avgConservation'] = data['conservationScore'].apply(
        lambda lst: np.mean(lst) if lst else np.nan
    )
    return data

def summarize_average_scores(data):
    # Drop any missing values (NaN) in the average scores.
    avg_scores = data['avgConservation'].dropna()
    mean_score = np.mean(avg_scores)
    median_score = np.median(avg_scores)
    std_score = np.std(avg_scores)
    return mean_score, median_score, std_score

def statistical_comparison(ion_scores, metabolite_scores):
    # Perform a t-test (assuming normality)
    t_stat, p_value_ttest = stats.ttest_ind(ion_scores, metabolite_scores, equal_var=False)
    
    # Perform a Mann-Whitney U test (non-parametric)
    u_stat, p_value_u = stats.mannwhitneyu(ion_scores, metabolite_scores)
    
    return t_stat, p_value_ttest, u_stat, p_value_u

# Step 5: Visualize the Results
def visualize_results(ion_scores, metabolite_scores):
    plt.figure(figsize=(10, 6))
    
    # Create a box plot to compare the conservation scores between ions and metabolites.
    # Note: Use 'tick_labels' instead of 'labels' if using Matplotlib 3.9+.
    plt.boxplot([ion_scores, metabolite_scores], labels=['Ions', 'Metabolites'])
    
    plt.title('Conservation Scores for Ions and Metabolites')
    plt.ylabel('Average Conservation Score')
    plt.show()
def visualize_results_violin(ion_scores, metabolite_scores):
    plt.figure(figsize=(10, 6))

    # Convert data into a format suitable for Seaborn
    data = {
        'Conservation Score': list(ion_scores) + list(metabolite_scores),
        'Group': ['Ions'] * len(ion_scores) + ['Metabolites'] * len(metabolite_scores)
    }
    df = pd.DataFrame(data)

    # Create a violin plot
    sns.violinplot(x='Group', y='Conservation Score', data=df, inner="quartile", palette="muted")

    plt.title('Distribution of Conservation Scores for Ions and Metabolites')
    plt.ylabel('Average Conservation Score')
    plt.xlabel('Ligand Type')
    plt.show()

def visualize_both_results(ion_scores, metabolite_scores):
    # Create a figure with 2 subplots (side by side)
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    
    # First plot: Boxplot for Ion vs. Metabolite scores
    axes[0].boxplot([ion_scores, metabolite_scores], labels=['Ions', 'Metabolites'])
    axes[0].set_title('Conservation Scores for Ions and Metabolites')
    axes[0].set_ylabel('Average Conservation Score')

    # Second plot: Violin plot for Ion vs. Metabolite scores
    data = {
        'Conservation Score': list(ion_scores) + list(metabolite_scores),
        'Group': ['Ions'] * len(ion_scores) + ['Metabolites'] * len(metabolite_scores)
    }
    df = pd.DataFrame(data)

    sns.violinplot(x='Group', y='Conservation Score', data=df, inner="quartile", palette="muted", ax=axes[1])
    axes[1].set_title('Distribution of Conservation Scores for Ions and Metabolites')
    axes[1].set_ylabel('Average Conservation Score')
    axes[1].set_xlabel('Ligand Type')

    # Adjust layout for better spacing
    plt.tight_layout()

    # Show both plots
    plt.show()

if __name__ == "__main__":
    file_path = FILE_PATHS.GET_ANALYSIS_RESULT_FILE("target_bacteria")  # Replace with your actual file path

    # Load and categorize the data.
    ion_data, metabolite_data = load_and_categorize_data(file_path)

    # Compute conservation scores for ions and metabolites.
    ion_data = compute_conservation_scores(ion_data)
    metabolite_data = compute_conservation_scores(metabolite_data)

    # Compute the average conservation score per entry.
    ion_data = compute_average_per_entry(ion_data)
    metabolite_data = compute_average_per_entry(metabolite_data)

    # Summarize conservation scores using the per-entry averages.
    ion_summary = summarize_average_scores(ion_data)
    metabolite_summary = summarize_average_scores(metabolite_data)

    # Print the summaries.
    print(f"Ion Conservation Score Summary: Mean = {ion_summary[0]}, Median = {ion_summary[1]}, Std Dev = {ion_summary[2]}")
    print(f"Metabolite Conservation Score Summary: Mean = {metabolite_summary[0]}, Median = {metabolite_summary[1]}, Std Dev = {metabolite_summary[2]}")

    # For statistical comparisons, extract the average scores.
    ion_avg_scores = ion_data['avgConservation'].dropna()
    metabolite_avg_scores = metabolite_data['avgConservation'].dropna()

    #check distribution
    shapiro_ion = stats.shapiro(ion_avg_scores)
    shapiro_metabolite = stats.shapiro(metabolite_avg_scores)

    print(f"Shapiro-Wilk Test for Ion Scores: W={shapiro_ion[0]}, p={shapiro_ion[1]}")
    print(f"Shapiro-Wilk Test for Metabolite Scores: W={shapiro_metabolite[0]}, p={shapiro_metabolite[1]}")

    # Perform statistical comparison.
    t_stat, p_value_ttest, u_stat, p_value_u = statistical_comparison(ion_avg_scores, metabolite_avg_scores)

    # Print statistical test results.
    print(f"T-test: t-statistic = {t_stat}, p-value = {p_value_ttest}")
    print(f"Mann-Whitney U test: U-statistic = {u_stat}, p-value = {p_value_u}")

    # Visualize the results.
    # visualize_results(ion_avg_scores, metabolite_avg_scores)
    # visualize_results_violin(ion_avg_scores, metabolite_avg_scores)
    visualize_both_results(ion_avg_scores, metabolite_avg_scores)
