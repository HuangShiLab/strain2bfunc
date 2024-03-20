###Strain2bFunc installer
###Updated at March 1, 2024
###Code by: ZHANG Yufeng
#!/bin/bash
##Users can change the default environment variables configuration file here
if [[ $SHELL = '/bin/zsh' ]];
then
        PATH_File=~/.zshrc
        if [ ! -f "$PATH_File" ]
        then
                PATH_File=~/.zsh_profile
                if [ ! -f "$PATH_File" ]
                then
                        touch $PATH_File
                fi
        fi
else
        PATH_File=~/.bashrc
        if [ ! -f "$PATH_File" ]
        then
                PATH_File=~/.bash_profile
                if [ ! -f "$PATH_File" ]
                then
                        touch $PATH_File
                fi
        fi

fi
Strain2bFunc=`pwd`
Sys_ver=`uname`
###Checking that environment variable of Strain2bFunc exists###
Check_old_strain2bfunc=`grep "export Strain2bFunc"  $PATH_File|awk -F '=' '{print $1}'`
Check_old_path=`grep "Strain2bFunc/bin"  $PATH_File |sed 's/\(.\).*/\1/' |awk '{if($1!="#"){print "Ture";}}'`
Add_Part="####DisabledbyStrain2bFunc####"
echo "**Strain2bFunc Installation**"
echo "**version 1.0**"

mkdir -p bin

###Build source code for src package###
if [ -f "Makefile" ]
   then
       echo -e "\n**Strain2bFunc src package**"
       make
       echo -e "\n**Build Complete**"
fi

chmod +x bin/Strain2bFunc-pipeline

###Configure environment variables###

if [ "$Check_old_strain2bfunc" != "" ]
   then
      Checking=`grep ^export\ Strain2bFunc  $PATH_File|awk -F '=' '{print $2}'`
      if [ "$Checking" != "$Strain2bFunc" ]
         then
         if [ "$Sys_ver" = "Darwin" ]
            then
            sed -i "" "s/^export\ Strain2bFunc/$Add_Part\ &/g" $PATH_File
            sed -i "" -e "`grep -n "$Add_Part" $PATH_File | cut -d ":" -f 1 | head -1` a\ 
            export\ Strain2bFunc=$Strain2bFunc
            " $PATH_File
         else
             sed -i "s/^export\ Strain2bFunc/$Add_Part\ &/g" $PATH_File
             sed -i "/$Add_Part\ export\ Strain2bFunc/a export\ Strain2bFunc=$Strain2bFunc" $PATH_File
         fi
     fi    
elif [ "$Check_old_strain2bfunc" = "" ]
    then
      echo "export Strain2bFunc="${Strain2bFunc} >> $PATH_File
fi

if [ "$Check_old_path" = "" ]
    then
      echo "export PATH=\$PATH:\$Strain2bFunc/bin" >> $PATH_File
fi

###Source the environment variable file###
source $PATH_File
echo -e "\n**Environment Variables Configuration Complete**"

echo -e "\n**Download database**"
cd database
perl Download_2bRADTagDB_GTDB.pl
cd ..

###End
echo -e "\n**Strain2bFunc Installation Complete**"
echo -e "\n**An example dataset with demo script is available in \"example\"**"
