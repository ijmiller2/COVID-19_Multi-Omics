
# Import create_engine
from sqlalchemy import create_engine, MetaData, Table, select

# Create an engine that connects to the census.sqlite file: engine
engine = create_engine('sqlite:///data/SQLite Database/Covid-19 Study DB.sqlite')

# Establish connection
connection = engine.connect()

# Print table names
print(engine.table_names())

# Instantiate metadata table
metadata = MetaData()

# Reflect the census table from the engine: census
metabolomics_measurements = Table('metabolomics_measurements', metadata, autoload=True, autoload_with=engine)

# Print the column names
print(metabolomics_measurements.columns.keys())

# Print full metadata of metabolomics_measurements
print(repr(metadata.tables['metabolomics_measurements']))

# SQL alchemy select all columns from metabolomics_measurements
stmt = select([metabolomics_measurements])
results = connection.execute(stmt).fetchall()

# Direct SQL query
stmt = 'SELECT * FROM metabolomics_measurements'
result_proxy = connection.execute(stmt)
results = result_proxy.fetchall()
