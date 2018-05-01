sleep 2s
echo "Thankyou for using ESMValTool for your diagnostic needs;"
echo "please pay attention to the following in-flight safety instructions"
sleep 2s
echo "========================================================="
echo "This is a script to run a quick ESMValTool test"
echo "========================================================="
echo "It will create an ESMValTool_QuickTest dir in the"
echo "$USER's $HOME directory; output will be stored there"
echo "so make sure you have enough space"
echo "It will get two files:"
echo " - config.yml : where the user's configuration is stored"
echo " - namelist.yml : where the run configuration is stored"
echo "If you get stuck email V at valeriu.predoi@ncas.ac.uk"
echo "Also remeber to check out the documentation"
echo "http://esmvaltool.readthedocs.io/en/refactoring_backend/"
echo "Have phun!"
echo "========================================================="
echo "Creating ESMValTool_QuickTest dir..."
mkdir -p ~/ESMValTool_QuickTest
cd ~/ESMValTool_QuickTest
echo "Getting config.yml from /home/users/valeriu/ESMValTool_KIT ..."
cp /home/users/valeriu/ESMValTool_KIT/config.yml .
echo "Getting namelist.yml from /home/users/valeriu/ESMValTool_KIT ..."
cp /home/users/valeriu/ESMValTool_KIT/namelist.yml .
echo "NOTE: data paths are set in config.yml;"
echo "      data lives in /home/users/valeriu/ESMValTool_KIT;"
echo "      please don't change that unless you have your own data!"
echo "OK let's run this thing..."
sleep 2s
echo "Running ESMValTool now..."
sleep 3s
export PATH=/home/users/valeriu/anaconda_users/envs/esmvaltool/lib/python3.6/site-packages:$PATH
export PATH=/home/users/valeriu/anaconda_users/envs/esmvaltool/bin:$PATH
/home/users/valeriu/anaconda_users/envs/esmvaltool/bin/esmvaltool -c config.yml -n namelist.yml
sleep 1s
cp /home/users/valeriu/ESMValTool_KIT/README.md .
sleep 1s
echo "Done...now read the README.md file that will explain what you have just done and you"
echo "will be able to understand the ESMValTool workflow process related to this quick test"
echo "Cheers! (btw you are still in the ESMValTool_QuickTest directory, don't panic :)"

