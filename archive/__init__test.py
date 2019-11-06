# import flask, render_template to render the html pages, and request to be able to handle post requests between html pages
from flask import Flask, render_template, request
import jinja2
import numpy as np


# import RandomForestClassifier to be able to evaluate models
from sklearn.ensemble import RandomForestClassifier

# import joblib to be able to load saved sklearn models
from sklearn.externals import joblib

# create a flask object
app = Flask(__name__)

htmlTitle = "Arindam Paul's"
template_text = "Testing"

@app.route('/')
def homepage():
	# this command returns the file and parameters that are rendered by the html file index.html
	return render_template("index.html", htmlTitle = htmlTitle, template_text = template_text)

if __name__ == "__main__":
	# to turn debug mode on use the following app.run(debug=True)
	app.run(debug=True)
