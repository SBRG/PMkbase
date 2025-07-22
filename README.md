
# Deploying the PMkbase/OmnilogDB Flask App on AWS EC2 (No SSL, IP Access Only)

This guide explains how to deploy the Flask-based PMkbase website on an Ubuntu EC2 instance using **Gunicorn** as the WSGI server and **Nginx** as the reverse proxy, without any SSL or domain configuration. The app will be accessible using your EC2 instance’s **public IP address**.

---

## 1. Prerequisites

1. An Ubuntu EC2 instance (tested with Ubuntu 20.04/22.04).  
2. SSH access with a user like `ubuntu`.  
3. Your Flask project directory uploaded to `/home/ubuntu/OmnilogDB`.  
4. A `requirements.txt` file in the project root.  
5. Security group configured to allow inbound traffic on ports `22 (SSH)`, `80 (HTTP)`.

---

## 2. Initial Server Setup

Update packages and install dependencies:

```bash
sudo apt update && sudo apt upgrade -y

# Install Python 3, pip, and virtualenv
sudo apt install python3 python3-pip python3-venv -y

# Install Nginx
sudo apt install nginx -y
```

---

## 3. Clone the repository and set Up the Python Virtual Environment

```bash

git clone https://github.com/SBRG/PMkbase.git
cd /home/ubuntu/PMkbase

# Create and activate the virtual environment
python3 -m venv pmkbase_env
source pmkbase_env/bin/activate

# Install dependencies
pip install -r requirements.txt
```

---

## 4. Test Run with Gunicorn

Before configuring systemd, test that Gunicorn can run your app:

```bash
# From the project root
gunicorn --workers 4 --bind 0.0.0.0:8000 app:app
```

- Replace `app:app` with `wsgi:app` if your app entry point is in `wsgi.py`.  
- Visit `http://<EC2-PUBLIC-IP>:8000` in your browser to verify it runs.  
- Stop it (`CTRL+C`) once verified.

---

## 5. Create a systemd Service for Gunicorn

Create the service file:

```bash
sudo nano /etc/systemd/system/pmkbase.service
```

Paste:

```ini
[Unit]
Description=Gunicorn instance for PMkbase
After=network.target

[Service]
User=ubuntu
Group=www-data
WorkingDirectory=/home/ubuntu/OmnilogDB
ExecStart=/home/ubuntu/myenv/bin/gunicorn -b localhost:8000 -w 20 --timeout 3600 app:app
Restart=always

[Install]
WantedBy=multi-user.target
```

Save and reload systemd:

```bash
sudo systemctl daemon-reload
sudo systemctl enable pmkbase
sudo systemctl start pmkbase

# Verify it's running
sudo systemctl status pmkbase
```

---

## 6. Configure Nginx

Edit the default Nginx configuration:

```bash
sudo nano /etc/nginx/sites-available/default
```

Replace contents with:

```nginx
upstream flaskpmkbase {
    server 127.0.0.1:8000;
}

server {
    listen 80 default_server;
    listen [::]:80 default_server;
    server_name _;

    root /var/www/html;

    index index.html index.htm index.nginx-debian.html;

    location / {
        proxy_pass http://flaskpmkbase;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
    }
}
```

Test and restart Nginx:

```bash
sudo nginx -t
sudo systemctl restart nginx
```

---

## 7. Managing the Service

Useful commands:

```bash
# Check Gunicorn service status
sudo systemctl status pmkbase

# Restart the service
sudo systemctl restart pmkbase

# Stop the service
sudo systemctl stop pmkbase

# Check if something is using port 8000
sudo lsof -i :8000

# Check Nginx status
sudo systemctl status nginx
```

To debug and run manually:

```bash
source myenv/bin/activate
gunicorn --workers 10 --timeout 3600 app:app
```

---

## 8. Access the App

Open your browser and go to:

```
http://<EC2-PUBLIC-IP>
```

Your Flask app should now be live.
