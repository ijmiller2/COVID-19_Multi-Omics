import matplotlib as mpl
mpl.use('Agg')
from optparse import OptionParser
import seaborn as sns
from matplotlib import pyplot as plt
import sys
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import numpy as np
from os.path import join
from sklearn.linear_model import ElasticNet, SGDRegressor
from sklearn.model_selection import GridSearchCV
from sklearn import metrics
import json

PARAMETERS = {
    'alpha': [
        0.001, 
        0.01, 
        0.1, 
        1.0,
        10.0,
        100.0
    ],
    'l1_ratio': [
        0.1, 
        0.2, 
        0.3, 
        0.4, 
        0.5, 
        0.6, 
        0.7, 
        0.8, 
        0.9, 
        1.0
    ]
}

MAX_ITER = 500000

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--subset", help="One of 'COVID', 'NONCOVID'")
    parser.add_option("-z", "--drop_zero_days", action='store_true', help='Drop patients with zero hospital-free days')
    parser.add_option("-l", "--log", action="store_true", help="Take log1 of data")
    parser.add_option("-o", "--output_prefix", help="Output file")
    (options, args) = parser.parse_args()

    expr_f = args[0]
    meta_f = args[1]
    prefix = options.output_prefix

    expr_df = pd.read_csv(expr_f, sep='\t', index_col=0)
    meta_df = pd.read_csv(meta_f, sep='\t')
    meta_df = meta_df.set_index('Albany_sampleID')
    print(meta_df)

    # Remove patients for which the 28 days has not elapsed
    no_expression_data = set(meta_df.index) - set(expr_df.columns)
    meta_df = meta_df.drop(no_expression_data)

    # Rmoeve non-COVID patients
    if options.subset == 'COVID':
        meta_df = meta_df.loc[meta_df['COVID'] == 1]
    elif options.subset == 'NONCOVID':
        meta_df = meta_df.loc[meta_df['COVID'] == 0]
       
    if options.drop_zero_days:
        meta_df = meta_df.loc[meta_df['Hospital_free_days_45'] > 0]
    print(meta_df)

    # Align metadata table and expression table
    expr_df = expr_df[meta_df.index]
    meta_df = meta_df.loc[expr_df.columns]

    hospital_free = np.array(meta_df['Hospital_free_days_45'])
    is_male = [
        int(g == 'M')
        for g in meta_df['Gender']
    ]
    print(is_male)
    is_female = [
        int(g == 'F')
        for g in meta_df['Gender']
    ]
    age = np.array(meta_df['Age_less_than_90'])
    X_clinical = np.array([
        is_male,
        is_female,
        age
    ]).T
   
    print(expr_df)

    X = np.array(expr_df)
    X = X.T
    if options.log:
        X = np.log(X+1)

    # Add clinical data
    print('Shape of gene expression matrix: {}'.format(X.shape))
    X = np.hstack([X, X_clinical])
    print('Shape after adding clinical data: {}'.format(X.shape))

    model = ElasticNet(max_iter=MAX_ITER)
    print('Performing grid-search cross-validation...')
    clf = GridSearchCV(
        model, 
        PARAMETERS, 
        scoring='neg_mean_squared_error'
    )
    clf.fit(X, hospital_free)
    print('done.')
    print(sorted(clf.cv_results_.keys()))
    print(clf.cv_results_['param_alpha'])
    print(clf.cv_results_['param_l1_ratio'])
    print(clf.cv_results_['mean_test_score'])

    # Select the best-performing parameters
    best_params = max(
        zip(
            clf.cv_results_['param_alpha'], 
            clf.cv_results_['param_l1_ratio'], 
            clf.cv_results_['mean_test_score']
        ),
        key=lambda x: x[2]
    )
    best_alpha = best_params[0]
    best_ratio = best_params[1]
    print('Max (alpha, ratio): ({}, {})'.format(best_alpha, best_ratio))

    # Fit the final model
    model = ElasticNet(
        alpha=best_alpha, 
        l1_ratio=best_ratio, 
        max_iter=MAX_ITER
    )
    model.fit(X, hospital_free)
    coeffs = model.coef_

    feats = list(expr_df.index) + ['is_male', 'is_female', 'age']
    print(feats)

    # Select positive coefficients
    up_model_coeffs = [
        (feat, coef)
        for feat, coef in zip(feats, coeffs)
        if coef > 0.0
    ]
    up_model_coeffs_df = pd.DataFrame(
        data=up_model_coeffs, 
        columns=['feature', 'coefficient']
    )
    up_model_coeffs_df = up_model_coeffs_df.set_index('feature')
    up_model_coeffs_df = up_model_coeffs_df.sort_values(
        by='coefficient', 
        axis=0, 
        ascending=False
    )

    # Write output file
    up_model_coeffs_df.to_csv(
        '{}.positive_features.tsv'.format(prefix), 
        sep='\t'
    )
    print("{} positive-zero features.".format(len(up_model_coeffs)))
    
    # Select negative coefficients
    down_model_coeffs = [
        (feat, coef)
        for feat, coef in zip(feats, coeffs)
        if coef < 0.0
    ]
    down_model_coeffs_df = pd.DataFrame(
        data=down_model_coeffs, 
        columns=['feature', 'coefficient']
    )
    down_model_coeffs_df = down_model_coeffs_df.set_index('feature')
    down_model_coeffs_df = down_model_coeffs_df.sort_values(
        by='coefficient', 
        axis=0, 
        ascending=False
    )

    # Write output file
    down_model_coeffs_df.to_csv(
        '{}.negative_features.tsv'.format(prefix), 
        sep='\t'
    )
    print("{} negative-zero features.".format(len(down_model_coeffs)))

    preds = model.predict(X)
    mae = metrics.mean_absolute_error(hospital_free, preds)
    with open('{}.model_fit.json'.format(prefix), 'w') as f:
        json.dump(
            {
                'Mean absolute error': mae
            },
            f,
            indent=True
        )


if __name__ == "__main__":
    main()
