import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

from app import app
from apps import pca, differential_expression, linear_regression, clustergrammer

app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content')
])

server = app.server

@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):

    if pathname == '/':
        return pca.layout

    elif pathname == '/pca':
        return pca.layout

    elif pathname == '/differential_expression':
        return differential_expression.layout

    elif pathname == '/linear_regression':
        return linear_regression.layout
    
    elif pathname == '/clustergrammer':
        return clustergrammer.layout

    else:
        return '404'

if __name__ == '__main__':
    app.run_server(
        debug=True,
        host='0.0.0.0')
