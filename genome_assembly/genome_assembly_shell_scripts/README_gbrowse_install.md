GBrowse-cloud-server
====================

Quickly deploy a GBrowse server in the cloud

**Step 1. Account setup**

The are numerous places to you can host your image, but I'll be using [DigitalOcean](https://www.digitalocean.com/) in this example. You'll need to go to that site and create an account, then select the type of instance to create. I chose `Ubuntu 12.04 x64` and the smallest plan. After you get an email with how to connect to your server, you will want to create a strong password by typing `passwd`. Then, you will need to start installing the dependencies.

**Step 2. Installing dependencies**

* Update your OS.

```bash
sudo apt-get update
```

* Then, upgrade.

```bash
sudo apt-get self-upgrade
```

* Install build tools like `gcc`, `make`, etc.

```bash
sudo apt-get install build-essential
```

* Install `git`.

```bash
sudo apt-get install git
```

**Step 3. Download GBrowse-cloud-server**
 
```bash
git clone https://github.com/sestaton/GBrowse-cloud-server
```

* Change to the project direcotry.

```bash
cd GBrowse-cloud-server
```

**Step 4. Run the install script to install BioPerl, GBrowse, and all dependencies**
```bash
bash gbrowse-cloud-server_install.sh
```

This final step is interactive and will take you through a number of questions. The results, including GBrowse configuration details,
will be logged for future reference.
