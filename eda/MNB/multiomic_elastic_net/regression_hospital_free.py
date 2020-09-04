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
from sklearn.linear_model import ElasticNet, Ridge
from sklearn.model_selection import GridSearchCV
from sklearn import metrics
from sklearn.preprocessing import scale
import json

PARAMETERS = {
    'elastic_net': {
        'alpha': [
            0.001, 
            0.01, 
            0.1, 
            1.0,
            10.0,
            100.0
        ],
        'l1_ratio': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    },
    'ridge': {
        'alpha': [
            0.001, 
            0.01, 
            0.1, 
            1.0,
            10.0,
            100.0
        ]
    }
}

MAX_ITER = 200000

META_COLS = [
    'sample_id',
    'Sample_label',
    'COVID',
    'Hospital_free_days',
    'Hospital_free_days_45',
    'Age_less_than_90',
    'Gender',
    'ICU_1',
    'APACHEII',
    'Charlson_score',
    'Mech_Ventilation',
    'Vent_free_days',
    'color_by',
    'DM',
    'Sample_label',	
    'Ferritin_NG/ML',	
    'CRP_MG/L',	
    'DDIMER_mg/L_FEU',	
    'Procalcitonin_NG/ML',
    'Lactate_MMol/L',	
    'Fibrinogen',	
    'SOFA'
]

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--subset", help="One of 'COVID', 'NONCOVID'")
    parser.add_option("-a", "--algorithm", help="Model: elastic_net, ridge")
    parser.add_option("-r", "--remove_unknowns", action="store_true", help="Remove unknown columns") 
    parser.add_option("-o", "--output_prefix", help="Output file")
    (options, args) = parser.parse_args()

    data_f = args[0]
    if options.algorithm:
        algo = options.algorithm
    else:
        algo = 'elastic_net'
    print('Running regression with {} penalty'.format(algo))
    prefix = options.output_prefix

    df = pd.read_csv(data_f, sep='\t')
    df = df.set_index('Albany_sampleID')

    # Select only COVID-patients
    df = df.loc[df['COVID'] == 1]
    meta_cols = set(META_COLS) & set(df.columns)
    meta_df = df[meta_cols]
    measurements_df = df.drop(meta_cols, axis=1)

    if options.remove_unknowns:
        keep_cols = [
            col
            for col in measurements_df.columns
            if 'unknown' not in col.lower()
        ]
        print('Keeping {}/{} columns after removing "unknonw"'.format(len(keep_cols), len(measurements_df.columns)))
        measurements_df = measurements_df[keep_cols]

    print(list(measurements_df.columns))
    print(measurements_df)
    hospital_free = np.array(meta_df['Hospital_free_days_45'])
    is_male = [
        int(g == 'M')
        for g in meta_df['Gender']
    ]
    is_female = [
        int(g == 'F')
        for g in meta_df['Gender']
    ]
    age = np.array(meta_df['Age_less_than_90'])
    
    # Build covariate matrix
    X_clinical = np.array([
        is_male,
        is_female,
        age
    ]).T
    X = np.array(measurements_df, dtype=np.float64)
    X = scale(X)
    print(X)

    # Add clinical data
    print('Shape of measurements matrix: {}'.format(X.shape))
    X = np.hstack([X, X_clinical])
    print('Shape after adding clinical data: {}'.format(X.shape))
    print(X)

    params = PARAMETERS[algo]
    if algo == 'elastic_net':
        model = ElasticNet(max_iter=MAX_ITER)
    elif algo == 'ridge':
        model = Ridge(max_iter=MAX_ITER)

    print('Performing grid-search cross-validation...')
    clf = GridSearchCV(model, params, scoring='neg_mean_squared_error')
    clf.fit(X, hospital_free)
    print('done.')
    print('Mean test scores:')
    print(clf.cv_results_['mean_test_score'])

    # Select best parameters
    if algo == 'elastic_net':
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
        model = ElasticNet(alpha=best_alpha, l1_ratio=best_ratio, max_iter=MAX_ITER)
    elif algo == 'ridge':
        best_params = max(
            zip(clf.cv_results_['param_alpha'], clf.cv_results_['mean_test_score']),
            key=lambda x: x[1]
        )
        best_alpha = best_params[0]
        print('Max alpha: {}'.format(best_alpha))
        model = Ridge(alpha=best_alpha)
    
    # Fit model on full dataset
    model.fit(X, hospital_free)
    coeffs = model.coef_

    # Feature names
    feats = list(measurements_df.columns) + ['is_male', 'is_female', 'age']

    # Select positive features
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
