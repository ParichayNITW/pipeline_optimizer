from flask import Flask, render_template, request, redirect, url_for
from solver import run_pipeline_optimization
import os

app = Flask(__name__)
app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', 'supersecret')

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        # Extract and validate form values
        try:
            data = {key: float(request.form[key]) for key in [
                'flow', 'kin_visc', 'density',
                'sfc_jamnagar', 'sfc_rajkot', 'sfc_surendranagar',
                'rate_dra', 'price_hsd'
            ]}
        except (ValueError, KeyError) as e:
            return render_template('form.html', error='Invalid input: please enter numeric values.', data=request.form)
        # Call the solver (exactly your code)
        try:
            results = run_pipeline_optimization(**data)
        except Exception as e:
            # Capture any solver or NEOS errors
            return render_template('form.html', error=f"Solver error: {e}", data=request.form)
        return render_template('results.html', results=results)
    return render_template('form.html')

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port)
