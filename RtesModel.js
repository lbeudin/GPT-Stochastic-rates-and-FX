// Model parameters (default values)
const models = {
    vasicek: { a: 0.1, b: 0.05, sigma: 0.02 },
    cir: { a: 0.1, b: 0.05, sigma: 0.02 },
    dothan: {  b: 0.05, sigma: 0.02 },
    bdt: {  b: 0.05, sigma: 0.02 },
    holee: { a:0.01, sigma: 0.02 },
    hullwhite: { a: 0.1, sigma: 0.02 }
};

// Initialize Chart
let ctx = document.getElementById('rateChart').getContext('2d');
let chart = new Chart(ctx, {
    type: 'line',
    data: { labels: [], datasets: [{ label: 'Asset Value', data: [] }] },
    options: { responsive: true }
});

// Function to simulate Vasicek model
function simulateVasicek(a, b, sigma, T = 1, steps = 100, r0 = 0.03) {
    let dt = T / steps;
    let r = r0; // Initial rate
    let values = [r];
    for (let i = 0; i < steps; i++) {
        let dw = Math.sqrt(dt) * gaussianRandom();
        r = r + a * (b - r) * dt + sigma * dw;
        values.push(r);
    }
    return values;
}
function simulateCir(a, b, sigma, T = 1, steps = 100, r0 = 0.03) {
    let dt = T / steps;
    let r = r0; // Initial rate
    let values = [r];
    for (let i = 0; i < steps; i++) {
        let dw = Math.sqrt(dt) * gaussianRandom();
        r = r + a * (b - r) * dt + sigma * Math.sqrt(r) * dw;
        values.push(r);
    }
    return values;
}

function simulateDothan(b, sigma, T = 1, steps = 100, r0 = 0.03) {
    let dt = T / steps;
    let r = r0; // Initial rate
    let values = [r];
    for (let i = 0; i < steps; i++) {
        let dw = Math.sqrt(dt) * gaussianRandom();
        r = r + b * r * dt + sigma * r * dw;
        values.push(r);
    }
    return values;
}
function simulateBdt(b, sigma, T = 1, steps = 100,r0 = 0.03) {
    let dt = T / steps;
    let r = r0; // Initial rate
    let values = [r];
    for (let i = 0; i < steps; i++) {
        let dw = Math.sqrt(dt) * gaussianRandom();
        r = r + b * r * dt + sigma  * dw;
        values.push(r);
    }
    return values;
}


function forwardRateCurve(t) {
    // Example: here ze should include the forward curve as function and will return the interpolated point for t
	// but I let chatGPT put its linar function
    return 0.03 + 0.002 * t;
}



// Function to compute the derivative of f(0, t) with respect to t
function forwardRateCurveDerivative(t, h=0.1) {
    return (forwardRateCurve(t+h)-forwardRateCurve(t))/h;
}

// Function to calculate b(t) based on the forward rate curve and its derivative
function b_t(t, a) {
    const f_0_t = forwardRateCurve(t); // f(0, t)
    const df_0_t_dt = forwardRateCurveDerivative(t); // Derivative of f(0, t)
    
    // b(t) = f(0, t) + (1 / a) * (df(0, t)/dt + a * f(0, t))
    return f_0_t + (1 / a) * (df_0_t_dt + a * f_0_t);
}

// Simulate the Hull-White model with time-dependent b(t)
function simulateHullWhite(a, sigma, T = 1, steps = 100, r0 = 0.03) {
    let dt = T / steps;
    let r = r0; // Initial rate
    let values = [r];

    for (let i = 0; i < steps; i++) {
        let t = i * dt;
        let b_t_value = b_t(t, a);  // Calculate b(t) at each step
        let dw = Math.sqrt(dt) * gaussianRandom();  // Brownian motion term
        r = r + a * (b_t_value - r) * dt + sigma * dw;  // Hull-White SDE
        values.push(r);  // Store simulated rate
    }

    return values;
}
// Simulate the Hull-White model with time-dependent b(t)
function simulateHoLee(a, sigma, T = 1, steps = 100, r0 = 0.03) {
    let dt = T / steps;
    let r = r0; // Initial rate
    let values = [r];

    for (let i = 0; i < steps; i++) {
        let t = i * dt;
        let b_t_value = b_t(t, a);  // Calculate b(t) at each step
        let dw = Math.sqrt(dt) * gaussianRandom();  // Brownian motion term
        r = r + (b_t_value ) * dt + sigma * dw;  
        values.push(r);  // Store simulated rate
    }

    return values;
}


// Gaussian Random Number generator
function gaussianRandom() {
    let u = Math.random(), v = Math.random();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}

// Event listener for model selection
document.getElementById('model').addEventListener('change', (event) => {
    const selectedModel = event.target.value;
    createSliders(selectedModel);
    updateChart(selectedModel);
});

// Function to create parameter sliders dynamically
function createSliders(model) {
    const slidersDiv = document.getElementById('sliders');
    slidersDiv.innerHTML = ''; // Clear existing sliders

    Object.keys(models[model]).forEach((param) => {
        const slider = document.createElement('input');
        slider.type = 'range';
        slider.min = '-5';
        slider.max = '5';
        slider.step = '0.01';
        slider.value = models[model][param];
        slider.id = param;

       const valueSpan = document.createElement('span');
        valueSpan.id = `${param.name}-value`;
        valueSpan.innerText = param.default; // Initialize with the default value

        const label = document.createElement('label');
        label.for = param;
        label.textContent = `${param}: `;

        slidersDiv.appendChild(label);
        slidersDiv.appendChild(slider);
        label.appendChild(valueSpan);  // Add value next to label

        // Add event listener to update chart on slider change
        slider.addEventListener('input', () => {
			            valueSpan.innerText = slider.value;
            models[model][param] = parseFloat(slider.value);
            updateChart(model);
        });
    });
}
// Update the chart based on selected model and parameters
function updateChart(model) {
    let values = [];
    if (model === 'vasicek') {
        const { a, b, sigma } = models[model];
        values = simulateVasicek(a, b, sigma);
    }
    // Handle other models similarly (CIR, Dothan, etc.)
    if (model === 'cir') {
        const { a, b, sigma } = models[model];
        values = simulateVasicek(a, b, sigma);
    }
	    if (model === 'dothan') {
        const {b, sigma } = models[model];
        values = simulateDothan(b, sigma);
    }
	if (model === 'bdt') {
        const {b, sigma } = models[model];
        values = simulateBdt(b, sigma);
    }
		if (model === 'hullwhite') {
        const {a, sigma } = models[model];
        values = simulateHullWhite(a, sigma);
    }
			if (model === 'holee') {
        const {a, sigma } = models[model];
        values = simulateHoLee(a, sigma);
    }
    // Update the chart
    chart.data.labels = Array.from({ length: values.length }, (_, i) => i);
    chart.data.datasets[0].data = values;
    chart.update();
}

// Initialize with Vasicek model
createSliders('vasicek');
updateChart('vasicek');
