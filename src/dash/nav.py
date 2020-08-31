
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html

COONLAB_LOGO="https://coonlabs.com/wp-content/uploads/2016/07/coon-logo-white.png"
navbar = dbc.NavbarSimple(

    [

    dbc.NavItem(dbc.NavLink(html.Span(
            "About",
            id="tooltip-about",
            style={"cursor": "pointer"}), href="https://www.medrxiv.org/content/10.1101/2020.07.17.20156513v1")),
        
    dbc.Tooltip(
        "View the Preprint",
        target="tooltip-about"
        ),

    dbc.DropdownMenu(
        children=[
            dbc.DropdownMenuItem("Meet the Team", header=True),
            dbc.DropdownMenuItem("UW–Madison", href="https://coonlabs.com/"),
            dbc.DropdownMenuItem("Albany Medical Center", href="https://www.amc.edu/Profiles/jaitova.cfm"),
            dbc.DropdownMenuItem("Morgridge Institute for Research", href="https://morgridge.org/research/regenerative-biology/bioinformatics/")
        ],
        nav=True,
        in_navbar=True,
        label="More",
    ),

    dbc.NavItem(html.Div(html.A(html.Img(src=COONLAB_LOGO,
            style={"height":"40px"}),
            href="https://coonlabs.com/"),
            className="d-none d-lg-block ml-4"))

    ],
    brand="NIH NATIONAL CENTER FOR QUANTITATIVE BIOLOGY OF COMPLEX SYSTEMS",
    brand_style={"font-size":"xx-large", "font-style":"italic"},
    brand_href="https://www.ncqbcs.com/",
    color="#5bc0de",
    dark=True,
    className="mt-1",
    fluid=True
)