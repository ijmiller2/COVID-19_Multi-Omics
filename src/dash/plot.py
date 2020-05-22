
import plotly.express as px
import pandas as pd
from sklearn.decomposition import PCA

color_dict = {
                "COVID_ICU":"#D53E4F",
                "COVID_NONICU":"#FDAE61",
                "NONCOVID_ICU":"#74ADD1",
                "NONCOVID_NONICU":"#66C2A5",
                "Col5":"#F46D43",
                "Col6":"#5AAE61",
                "Col7":"#8073AC",
                "Col8":"#DE77AE",
                "Col9":"#9E0142",
                "Col10":"#F4A582",
                "Col11":"#2A4023",
                "Col12":"#2C0379"
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

def biomolecule_bar(combined_df, x, y, biomolecule_name):

    color_list = get_color_list(combined_df)

    fig = px.bar(combined_df, x=x, y=y, color=color_list, color_discrete_map=color_dict)

    fig.update_layout(
        title="Biomolecule Bar Plot",
        xaxis_title='Sample ID',
        yaxis_title='{} log2(LFQ)'.format(biomolecule_name),
        font=dict(
            family="Helvetica",
            size=18,
            color="#7f7f7f"))

    return fig

def boxplot(combined_df, biomolecule_name):

    color_list = get_color_list(combined_df)

    df = pd.DataFrame({'y':combined_df[biomolecule_name], 'color':color_list})

    fig = px.box(df, y="y", color="color", color_discrete_map=color_dict,
                points='all')
    fig.update_traces(quartilemethod="exclusive") # or "inclusive", or "linear" by default
    fig.update_layout(
        title="Biomolecule Box Plot",
        xaxis_title='Group',
        yaxis_title='{} log2(LFQ)'.format(biomolecule_name),
        font=dict(
            family="Helvetica",
            size=18,
            color="#7f7f7f"))


    return fig

def pca_plot(quant_df, combined_df):

    from sklearn.decomposition import PCA

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
                    size=[10]*df.shape[0],
                    color=color_list,
                    color_discrete_map=color_dict)
                    #symbol=shape_list,
                    #symbol_map={"COVID":'hexagram', "NONCOVID":'circle', "NA":'cross'})

    fig.update_layout(
        title="GC/MS Metabolomics PCA",
        xaxis_title='PC1 ({}%)'.format(round(100*pca.explained_variance_ratio_[0],1)),
        yaxis_title='PC2 ({}%)'.format(round(100*pca.explained_variance_ratio_[1],1)),
        font=dict(
            family="Helvetica",
            size=18,
            color="#7f7f7f")
        )

    return fig