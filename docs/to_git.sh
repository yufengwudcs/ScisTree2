#!/bin/bash

cd ../scistree2doc
cp -r ../docs .
git add .
git commit -m "update docs"
git push origin main
