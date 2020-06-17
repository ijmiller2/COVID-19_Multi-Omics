
import plotly.express as px
import pandas as pd
import numpy as np
from scipy.spatial import distance
from plotly import graph_objects
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import summary_table

color_dict = {
                "COVID_ICU":"#D53E4F",
                "COVID_NONICU":"#FDAE61",
                "NONCOVID_ICU":"#74ADD1",
                "NONCOVID_NONICU":"#66C2A5",
                "Male":"#F46D43",
                "Female":"#5AAE61",
                "Col7":"#8073AC",
                "Col8":"#DE77AE",
                "proteomics":"#9E0142",
                "lipidomics":"#F4A582",
                "metabolomics":"#2A4023",
                "transcriptomics":"#2C0379",
                "selected_biomolecule":"black"
                }

def get_color_list(combined_df):
    # from combined_df

    # get colors
    color_list = []
    for sample_id, row in combined_df.iterrows():

        ICU_1 = row['ICU_1']
        COVID = row['COVID']

        if pd.isnull(ICU_1):
            color = "Col12"

        elif ICU_1 == 1 and COVID == 1:
            color = color_dict["COVID_ICU"]
            color = "COVID_ICU"

        elif ICU_1 == 1 and COVID == 0:
            color = color_dict["NONCOVID_ICU"]
            color = "NONCOVID_ICU"

        elif ICU_1 == 0 and COVID == 1:
            color = color_dict["COVID_NONICU"]
            color = 'COVID_NONICU'

        elif ICU_1 == 0 and COVID == 0:
            color = color_dict["NONCOVID_NONICU"]
            color = "NONCOVID_NONICU"

        color_list.append(color)

    return color_list

def biomolecule_bar(combined_df, biomolecule_id, biomolecule_names_dict):

    biomolecule_name = biomolecule_names_dict[biomolecule_id]

    # sort the samples by group
    color_list = get_color_list(combined_df)
    combined_df['color_by'] = color_list
    combined_df['sample'] = combined_df.index
    combined_df.sort_values(by=['color_by', 'sample'], inplace=True)

    fig = px.bar(combined_df, x=[i for i in range(combined_df.shape[0])],
        y=combined_df[biomolecule_id],
        color=combined_df['color_by'],
        hover_data=['sample'],
        color_discrete_map=color_dict)

    fig.update_layout(
        title="{}".format(biomolecule_name),
        legend_title_text='Group',
        xaxis_title='Sample',
        yaxis_title='log2(LFQ) Value',
        font=dict(
            family="Helvetica",
            size=18,
            color="#7f7f7f"))

    return fig

def boxplot(combined_df, biomolecule_id, biomolecule_names_dict):

    biomolecule_name = biomolecule_names_dict[biomolecule_id]

    color_list = get_color_list(combined_df)

    df = pd.DataFrame({'y':combined_df[biomolecule_id], 'color':color_list,
                        'sample':combined_df.index})

    fig = px.box(df, y="y", color="color", color_discrete_map=color_dict,
                points='all', hover_data=['sample'])

    fig.update_traces(quartilemethod="exclusive") # or "inclusive", or "linear" by default
    fig.update_layout(
        title="{}".format(biomolecule_name),
        legend_title_text='Group',
        xaxis_title='Group',
        yaxis_title='log2(LFQ) Value',
        showlegend=False,
        font=dict(
            family="Helvetica",
            size=18,
            color="#7f7f7f"))


    return fig

def pca_scores_plot(combined_df, quant_value_range):

    from sklearn.decomposition import PCA

    quant_columns = combined_df.columns[:quant_value_range]
    quant_df = combined_df[quant_columns]

    pca = PCA(n_components = 10)
    PCA = pca.fit_transform(quant_df)

    PC1s = []
    PC2s = []
    for PCs in PCA:
        PC1 = PCs[0]
        PC2 = PCs[1]
        PC1s.append(PC1)
        PC2s.append(PC2)

    color_list = get_color_list(combined_df)

    df = pd.DataFrame({'x':PC1s, 'y':PC2s,
        'sample_id':combined_df.index.tolist(),
        'COVID':combined_df['COVID']})

    fig = px.scatter(df, x="x", y="y", hover_data=['sample_id'],
                    color=color_list,
                    color_discrete_map=color_dict,
                    symbol='COVID')

    fig.update_traces(marker=dict(size=15, opacity=0.8))

    fig.update_layout(
        title="Samples (n={})".format(quant_df.shape[0]),
        legend_title_text='Group',
        xaxis_title='PC1 ({}%)'.format(round(100*pca.explained_variance_ratio_[0],1)),
        yaxis_title='PC2 ({}%)'.format(round(100*pca.explained_variance_ratio_[1],1)),
        showlegend=False,
        font=dict(
            family="Helvetica",
            size=18,
            color="#7f7f7f")
        )

    return fig

def pca_loadings_plot(combined_df, quant_value_range, dataset_id, biomolecule_names_dict, ome_type_list):

    from sklearn.decomposition import PCA

    quant_columns = combined_df.columns[:quant_value_range]
    # # NOTE: For some reason, quant_df here ends up with larger shape...
    quant_df = combined_df[quant_columns]

    #  NOTE: There appear to be duplicate lipid names
    # All lipid features currently set to keep=1
    quant_df = quant_df.loc[:,~quant_df.columns.duplicated()]

    pca = PCA(n_components = 10)
    PCA = pca.fit_transform(quant_df)

    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)

    PC1_index = 0
    PC1_loadings = [x[PC1_index] for x in loadings]
    PC2_index = 1
    PC2_loadings = [y[PC2_index] for y in loadings]

    df = pd.DataFrame({'x':PC1_loadings, 'y':PC2_loadings,
        'biomolecule_id':quant_df.columns.tolist(),
        'standardized_name':[biomolecule_names_dict[i] for i in quant_df.columns.tolist()],
        'ome_type':ome_type_list})

    # downsample larger plots
    if df.shape[0] > 1000:
        keep_list = downsample_scatter_data_by_variance(quant_df)

        df_drop_list = []
        for index,row in df.iterrows():
            if not row['biomolecule_id'] in keep_list:
                df_drop_list.append(index)
        df = df.drop(df_drop_list)

    fig = px.scatter(df, x="x", y="y",
        hover_data=['biomolecule_id', 'standardized_name'],
        color="ome_type",
        color_discrete_map=color_dict)

    fig.update_traces(marker=dict(size=10, opacity=0.5))

    fig.update_layout(
        title="{} (n={})".format(dataset_id, quant_df.shape[1]),
        legend_title_text='Group',
        xaxis_title='Loadings on PC1',
        yaxis_title='Loadings on PC2',
        font=dict(
            family="Helvetica",
            size=18,
            color="#7f7f7f")
        )

    # only show the color legend with combined datasets
    #if not dataset_id=="Combined":
    #    fig.update_layout(showlegend=False)

    return fig

def downsample_scatter_data(df):

    # df should have, x, y, and ome_type columns

    # filter top n biomolecules on loadings plot, by distance from origin
    origin = (0,0)
    distance_list = []
    for index, row in df.iterrows():
        x = row['x']
        y = row['y']
        coordinates = (x, y)
        d = distance.euclidean(origin, coordinates)
        distance_list.append(d)

    df['distance_from_origin'] = distance_list

    distance_std = np.std(distance_list)
    downsample_range = distance_std

    drop_index_list = []
    for ome_type in list(set(df['ome_type'])):
        # drop 20 % of measurements for each ome, randomly subsample 50% of those
        #drop_row_num = round(df[df['ome_type'] == ome_type].shape[0] * 0.2)
        #drop_indices = df[df['ome_type'] == ome_type].\
        #    sort_values(by='distance_from_origin').\
        #    iloc[:drop_row_num].sample(frac=0.5, random_state=1).index.tolist()

        # randomly downsample half of data within one standard deviation from origin
        ome_df = df[(df['ome_type'] == ome_type) & (df['distance_from_origin'] < downsample_range)]
        drop_indices = ome_df.sample(random_state=1, frac=0.25).index.tolist()
        drop_index_list.extend(drop_indices)

    df = df.drop(drop_index_list)

    return df

def downsample_scatter_data_by_variance(df):
    # return top 1000 features by variance
    keep_list = df.std(axis=0).sort_values(ascending=False)[:1000]

    return keep_list

def downsample_volcano_data(df):

    keep_index_list = []
    for ome_type in list(set(df['ome_type'])):
        # keep data for each ome with top 1000 features by variance (std)
        ome_df = df[df['ome_type'] == ome_type].sort_values(by='std', ascending=False)
        #ome_df = df[df['ome_type'] == ome_type].sort_values(by='q_value', ascending=True)
        keep_indices = ome_df.index.tolist()[:2000]
        keep_index_list.extend(keep_indices)

    df = df.loc[keep_index_list]

    return df

def volcano_plot(volcano_df):

    #volcano_df.dropna(inplace=True)
    volcano_df = volcano_df.dropna()

    df = pd.DataFrame({'x':volcano_df['log2_FC'],
        'y':volcano_df['neg_log10_p_value'],
        'biomolecule_id':volcano_df['biomolecule_id'],
        'standardized_name':volcano_df['standardized_name'],
        'ome_type':volcano_df['ome_type'],
        'p_value':volcano_df['p_value'],
        'q_value':volcano_df['q_value'],
        'std':volcano_df['std']})

    df = downsample_volcano_data(df)

    fig = px.scatter(df, x="x", y="y",
    hover_data=['biomolecule_id', 'standardized_name', 'p_value', 'q_value'],
    opacity=0.5,
    size='y',
    color='ome_type',
    color_discrete_map=color_dict)

    #fig.update_traces(marker=dict(size=10, opacity=0.5))

    #confounders = ", ".join(volcano_df.iloc[0]['confounders'].split(";"))
    title = "COVID vs NONCOVID"

    fig.update_layout(
        title="{} (n={})".format(title, volcano_df.shape[0]),
        legend_title_text='Dataset',
        xaxis_title='Effect Size (log2 FC)',
        yaxis_title='Significance (-log10(Corrected P-value))',
        font=dict(
            family="Helvetica",
            size=18,
            color="#7f7f7f")
        )

    return fig

def correlation_scatter(combined_df, biomolecule_id, selected_groups,
    biomolecule_name, clinical_measurement):

    # shorten biomolecule name
    if len(biomolecule_name) > 15:
        biomolecule_name = biomolecule_name[:15] + ".."

    ## NOTE: See https://pythonplot.com/ for confidence interval example
    if clinical_measurement == "Gender":
        combined_df.replace("M", 0, inplace=True)
        combined_df.replace("F", 1, inplace=True)

    # drop samples with missing values for clinical measurement
    combined_df.replace('', np.nan, inplace=True)
    combined_df = combined_df.dropna(subset=[clinical_measurement, 'COVID'])

    color_list = get_color_list(combined_df)
    combined_df['group'] = color_list

    # drop by group
    for group in ['COVID_ICU', 'COVID_NONICU', 'NONCOVID_ICU', "NONCOVID_NONICU"]:
        if not group in selected_groups:
            combined_df = combined_df[combined_df['group'] != group]

    # set target and explanatory variables
    y_var = biomolecule_id
    x_var = clinical_measurement

    ### Run regression with statsmodels ###
    x=combined_df[x_var].astype(float)
    y=combined_df[y_var].astype(float)
    X = sm.add_constant(x)
    res = sm.OLS(y, X).fit()
    rsquared = round(res.rsquared, 3)
    p_val = '%.3E' % res.f_pvalue

    # get regression data for range of HFD values
    x_min = round(min(x),1)
    x_max = round(max(x), 1)
    x_range = np.arange(x_min,x_max + 0.1 ,0.1)
    #x_range = np.array([i for i in range(min(combined_df[biomolecule_id]), max(combined_df[biomolecule_id], 0.1))])
    X = sm.add_constant(x_range)
    out_of_sample_predictions = res.get_prediction(X)
    preds = out_of_sample_predictions.summary_frame(alpha=0.05)

    ###

    df = pd.DataFrame({'x':x,
        'y':y,
        'sample_id':combined_df.index.tolist(),
        'COVID':combined_df['COVID'],
        'group':combined_df['group']})

    fig = px.scatter(df, x="x", y="y", hover_data=['sample_id'],
                    color='group',
                    color_discrete_map=color_dict)

    fig.update_traces(marker=dict(size=15, opacity=0.8))

    # add regression data
    line_of_best_fit = graph_objects.Scatter({
    'mode' : 'lines',
    'x' : x_range,
    'y' : preds['mean'],
    'name' : 'Trend',
    'opacity' : 0.6,
    'line' : {
        'color' : 'black'
    }
    })

    #Add a lower bound for the confidence interval, white
    mean_ci_lower = graph_objects.Scatter({
        'mode' : 'lines',
        'x' : x_range,
        'y' : preds['mean_ci_lower'],
        'name' : 'Lower 95% CI',
        'showlegend' : False,
        'line' : {
            'color' : 'white'
        }
    })
    # Upper bound for the confidence band, transparent but with fill
    mean_ci_upper = graph_objects.Scatter( {
        'type' : 'scatter',
        'mode' : 'lines',
        'x' : x_range,
        'y' : preds['mean_ci_upper'],
        'name' : '95% CI',
        'fill' : 'tonexty',
        'line' : {
            'color' : 'white'
        },
        'fillcolor' : 'rgba(255, 127, 14, 0.3)'
    })

    fig.add_trace(line_of_best_fit)
    fig.add_trace(mean_ci_lower)
    fig.add_trace(mean_ci_upper)

    if len(x_var) > 10:
        x_var = x_var[:10] + ".."
    formula = "Biomolecule {} ~ {}".format(biomolecule_id, x_var)
    plot_title = "{}, R2: {}, p value: {}, n={}".format(formula, rsquared, p_val, combined_df.shape[0])

    fig.update_layout(
        title=plot_title,
        legend_title_text='Group',
        xaxis_title='{}'.format(clinical_measurement),
        yaxis_title='{} \nlog2 intensity'.format(biomolecule_name),
        showlegend=True,
        font=dict(
            family="Helvetica",
            size=18,
            color="#7f7f7f")
        )

    return fig
