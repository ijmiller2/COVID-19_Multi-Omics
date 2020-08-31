
## Notes on launching server

1. Transfer SQLITE database to 'data/SQLite Database/DATE'

2. Sync git repo
- `git checkout --force # to avoid conflicts from local changes`
- `git pull`

3. Update `index.py`
- `app.run_server(debug=False, host='0.0.0.0', port='8080')`
- Running on: 209.188.7.206:8080

4. Run with anaconda3:
- start screen: `screen`
- Or reattach to screen: `screen -ls` + `scree -r $screen_id`
- `/root/anaconda3/bin/python index.py`
- Detach from screen: `screen cntrl+a cntrl+d`

5. Running with gunicorn:
- `/root/anaconda3/bin/gunicorn index:server -b 0.0.0.0:8080 --timeout=180`

6. Running via https:
- `/root/anaconda3/bin/gunicorn index:server -b 0.0.0.0:8080 --timeout=180 --certfile covid-omics.app.crt --keyfile covid-omics.key`