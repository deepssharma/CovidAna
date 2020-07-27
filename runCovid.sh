#!/bin/bash

cd /Users/deepalisharma/git/COVID-19;
git pull;
cd /Users/deepalisharma/git/CovidAna;
ls -1 /Users/deepalisharma/git/COVID-19/csse_covid_19_data/csse_covid_19_daily_reports_us/*.csv > data.list;
echo "Results will be shown for the states: $2, $3, $4";
check=1;
if [ $check == $1 ]
then
    root -l -q make_covid_tree.C\("$1",\"$2\",\"$3\",\"$4\"\)
    root -l -q make_covid_tree.C\(0,\"$2\",\"$3\",\"$4\"\)
else
    root -l -q make_covid_tree.C\(0,\"$2\",\"$3\",\"$4\"\);
fi
