
import plotly.express as px
import pandas as pd
import numpy as np

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

        elif ICU_1 == "1" and COVID == "1":
            color = color_dict["COVID_ICU"]
            color = "COVID_ICU"

        elif ICU_1 == "1" and COVID == "0":
            color = color_dict["NONCOVID_ICU"]
            color = "NONCOVID_ICU"

        elif ICU_1 == "0" and COVID == "1":
            color = color_dict["COVID_NONICU"]
            color = 'COVID_NONICU'

        elif ICU_1 == "0" and COVID == "0":
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

    df = pd.DataFrame({'x':PC1s, 'y':PC2s, 'sample_id':combined_df.index.tolist()})

    fig = px.scatter(df, x="x", y="y", hover_data=['sample_id'],
                    color=color_list,
                    color_discrete_map=color_dict)

    fig.update_traces(marker=dict(size=20, opacity=0.8))

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
        'standardized_name':[biomolecule_names_dict[i] for i in quant_df.columns.tolist()]})

    fig = px.scatter(df, x="x", y="y",
        hover_data=['biomolecule_id', 'standardized_name'],
        color=ome_type_list,
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
