        // Store the parameters and formula for each model
        const models = {
            dupire: {
                formula: 'dSt = sigma(St, t) * St * dWt',
                description: 'The Dupire local volatility model assumes that the volatility is a deterministic function of both time and the asset price.',
                parameters: [
                    { name: 'σ', min: 0.1, max: 1, step: 0.01, default: 0.2 }
                ],
                simulate: (params) => simulateDupire(params)
            },
            ecir: {
                formula: 'dSt = (rd - rf - lambda(t) * vt) * St * dt + sqrt{vt} * St * dWt',
                description: 'The Extended CIR (ECIR) model extends the CIR model by allowing time-dependence in parameters, which is particularly suited for FX markets.',
                parameters: [
                    { name: 'κ  Mean reversion rate of Var', min: 0.01, max: 0.5, step: 0.01, default: 0.1 },
                    { name: 'θ Long-term Var', min: 0.01, max: 1, step: 0.01, default: 0.2 },
                    { name: 'σ Vol Var', min: 0.01, max: 1, step: 0.01, default: 0.2 }
                ],
                simulate: (params) => simulateECIR(params)
            },

			bates: {
                formula: 'dS_t = ',
                description: '',
                parameters: [
                    { name: 'κ Mean reversion rate of Var ', min: -0.5, max: 0.5, step: 0.01, default: 0.1 },
                    { name: 'θ Long-term Var', min: 0.01, max: 1, step: 0.01, default: 0.2 },
                    { name: 'σ Vol Var', min: 0.01, max: 10, step: 0.01, default: 0.2 },
                    { name: 'ρ Correl Asset and Var', min: 0.01, max: 1, step: 0.01, default: 0.2 },
                    { name: 'v0 init Var', min: 0.01, max: 1, step: 0.01, default: 0.3 },
                    { name: 'λ jump intensity avgJump/unit time', min: 0.0001, max: 1, step: 0.01, default: 0.05 },
                    { name: 'ν mean jump size distribution', min: 0.00001, max: 0.0010, step: 0.0001, default: 0.0002 },
                    { name: 'δ std deviation jump size distri', min: 0.000001, max: 1, step: 0.0001, default: 0.0002 }
                ],
                simulate: (params) => simulateBates(params)
            },
			
			merton: {
				formula: 'dSt = μ * St *dt + σ*St * dWt + (J - 1) St dNt',
				description: 'The Merton Jump Diffusion model incorporates sudden jumps into the Black-Scholes model, allowing for sudden shifts in FX rates. Jumps occur randomly with a log-normal distribution.',
				parameters: [
					{ name: 'μ Drift of the diffusion process', min: -0.5, max: 0.5, step: 0.01, default: 0.1 },
					{ name: 'σ  Volatility of the diffusion process', min: 0.01, max: 1, step: 0.01, default: 0.2 },
					{ name: 'λ Jump intensity (average number of jumps per unit time)', min: 0.01, max: 1, step: 0.01, default: 0.1 },
					{ name: 'ν Mean of the jump size in log-space (log-normal distribution)', min: -0.5, max: 0.5, step: 0.01, default: 0.1 },
					{ name: 'δ Standard deviation of the jump size in log-space', min: 0.01, max: 0.5, step: 0.01, default: 0.1 }
				],
				simulate: (params) => simulateMertonJump(params)
			},

			heston: {
				formula: 'dSt = μ * St * dt + sqrt{vt} St + dWt^S with dvt = κ (θ - vt) * dt + σ *sqrt{vt} * dWt^v',
				description: 'The Heston Stochastic Volatility Model allows volatility to vary over time with mean-reversion, capturing realistic market behaviors such as volatility clustering.',
				parameters: [
					{ name: 'μ Drift of the asset price', min: -0.5, max: 0.5, step: 0.01, default: 0.05 },
					{ name: 'κ Mean reversion rate of the variance', min: 0.01, max: 5, step: 0.01, default: 0.5 },
					{ name: 'θ Long-term mean of the variance', min: 0.01, max: 1, step: 0.01, default: 0.2 },
					{ name: 'σ Volatility of the variance (vol of vol).', min: 0.01, max: 1, step: 0.01, default: 0.3 },
					{ name: 'ρ  Correlation between asset price and variance', min: -1, max: 1, step: 0.01, default: 0.5 },
					{ name: 'v0 Initial variance', min: 0.01, max: 1, step: 0.01, default: 0.2 }
				],
				simulate: (params) => simulateHeston(params)
			},

			garch: {
				formula: 'St = S{t-1} exp(σt  Z_t)   σt^2 = ω + α * epsilon{t-1}^2 + β σ{t-1}^2',
				description: 'The GARCH(1,1) Model captures time-varying volatility by adjusting based on past returns and volatility, useful for modeling volatility clustering.',
				parameters: [
					{ name: 'ω vol level', min: 0.0001, max: 1, step: 0.0001, default: 0.01 },
					{ name: 'α sensitivity to recent returns', min: 0, max: 1, step: 0.01, default: 0.1 },
					{ name: 'β volatility persistence', min: 0, max: 1, step: 0.01, default: 0.8 }
				],
				simulate: (params) => simulateGARCH(params)
			},

			vg: {
				formula: 'St = S0 exp((μ - 0.5 σ^2 - θ) t + σ W{μ}(t) + θ G{ν}(t))',
				description: 'The Variance Gamma Model captures asset price behavior with heavy tails and skewness, incorporating jumps in a Gamma-distributed process.',
				parameters: [
					{ name: 'μ drift', min: -0.5, max: 0.5, step: 0.01, default: 0.05 },
					{ name: 'σ vol of the diffusion', min: 0.01, max: 1, step: 0.01, default: 0.2 },
					{ name: 'θ skewness', min: -1, max: 1, step: 0.01, default: 0 },
					{ name: 'ν Variance rate of the Gamma process, affecting the rate of jumps.', min: 0.01, max: 1, step: 0.01, default: 0.1 }
				],
				simulate: (params) => simulateVarianceGamma(params)
			},
			sabr: {
				formula: 'dSt = σt St^β  dWt^S d  &&  σt = α sigma_t  dW_t^σ',
				description: 'The SABR model captures stochastic volatility with a correlation between the asset price and volatility, often used to model FX options and other derivatives.',
				parameters: [
					{ name: 'β  elasticity of the volatility', min: 0, max: 1, step: 0.01, default: 0.5 },
					{ name: 'α vol of vol', min: 0.01, max: 1, step: 0.01, default: 0.3 },
					{ name: 'ρ correlation between price and volatility shocks', min: -1, max: 1, step: 0.01, default: 0.5 },
					{ name: 'σ0 init vol', min: 0.01, max: 1, step: 0.01, default: 0.2 }
				],
				simulate: (params) => simulateSABR(params)
			}


            // Add more models here with corresponding parameters
            // GARCH, SABR, 
        };

        // Function to update the formula, description, and sliders dynamically based on the selected model
        function updateModel() {
            const modelKey = document.getElementById('model-select').value;
            const model = models[modelKey];

            // Update formula
            document.getElementById('formula-content').innerText = model.formula;

            // Update description
            document.getElementById('description-content').innerText = model.description;

            // Create sliders for the parameters
            createSliders(model);
        }



        // Function to create sliders for each parameter
        function createSliders(model) {
            const sliderContainer = document.getElementById('slider-container');
            sliderContainer.innerHTML = ''; // Clear previous sliders

            model.parameters.forEach(param => {
                const div = document.createElement('div');
                div.className = 'slider-container';

                const label = document.createElement('label');
                label.innerText = `${param.name}: `;

                const valueSpan = document.createElement('span');
                valueSpan.id = `${param.name}-value`;
                valueSpan.innerText = param.default;

                const slider = document.createElement('input');
                slider.type = 'range';
                slider.min = param.min;
                slider.max = param.max;
                slider.step = param.step;
                slider.value = param.default;
                slider.id = `${param.name}-slider`;

                // Update the displayed value when the slider is adjusted
                slider.addEventListener('input', () => {
                    valueSpan.innerText = slider.value;
                    param.value = parseFloat(slider.value);
                    updateChart(model); // Update chart with new parameter values
                });

                label.appendChild(valueSpan);  // Add value next to label
                div.appendChild(label);
                div.appendChild(slider);
                sliderContainer.appendChild(div);
            });

            updateChart(model); // Initial chart with default values
        }

        // Function to update the chart dynamically based on model and parameters
        function updateChart(model) {
            const ctx = document.getElementById('chart').getContext('2d');

            // Simulate the FX model with current parameters
            const values = model.simulate(model.parameters);

            // Update the chart
            if (window.myChart) {
                window.myChart.destroy();
            }
            window.myChart = new Chart(ctx, {
                type: 'line',
                data: {
                    labels: Array.from({ length: 100 }, (_, i) => i), // 100 time steps
                    datasets: [{
                        label: 'Simulated FX Rate',
                        data: values,
                        borderColor: 'rgba(75, 192, 192, 1)',
                        borderWidth: 2,
                        fill: false
                    }]
                },
                options: {
                    responsive: true,
                    scales: {
                        x: {
                            title: { display: true, text: 'Time' }
                        },
                        y: {
                            title: { display: true, text: 'FX Rate' }
                        }
                    }
                }
            });
        }
		function simulateSABR(params) {
			const [beta, alpha, rho, sigma0] = params.map(p => p.value);
			let S = 1; // Initial FX rate
			let sigma = sigma0; // Initial volatility
			const dt = 0.01; // Time step
			const steps = 100; // Number of time steps
			const data = [S];

			for (let i = 0; i < steps; i++) {
				// Generate correlated Brownian increments
				const dW_S = Math.sqrt(dt) * gaussianRandom();
				const dW_sigma = Math.sqrt(dt) * gaussianRandom();
				const dW_corr = rho * dW_S + Math.sqrt(1 - rho ** 2) * dW_sigma; // Correlated Brownian motion

				// Volatility process
				const d_sigma = alpha * sigma * dW_corr;
				sigma = Math.max(sigma + d_sigma, 0); // Ensure non-negative volatility

				// Asset price process
				const dS = sigma * S ** beta * dW_S;
				S = S + dS;

				data.push(S);
			}

			return data;
		}

		function simulateHeston(params) {
			const [mu, kappa, theta, sigma, rho, v0] = params.map(p => p.value);
			let S = 1; // Initial FX rate
			let v = v0; // Initial variance
			const dt = 0.01; // Time step
			const steps = 100; // Number of time steps
			const data = [S];

			for (let i = 0; i < steps; i++) {
				// Generate correlated Brownian increments
				const dW_S = Math.sqrt(dt) * gaussianRandom();
				const dW_v = Math.sqrt(dt) * gaussianRandom();
				const dW_corr = rho * dW_S + Math.sqrt(1 - rho ** 2) * dW_v; // Correlated Brownian motion

				// Variance process (mean-reverting)
				const dv = kappa * (theta - v) * dt + sigma * Math.sqrt(Math.max(v, 0)) * dW_corr;
				v = Math.max(v + dv, 0); // Ensure variance stays non-negative

				// Asset price process
				const dS = mu * S * dt + Math.sqrt(v) * S * dW_S;
				S = S + dS;
				
				data.push(S);
			}

			return data;
		}
		
		function simulateGARCH(params) {
			const [omega, alpha, beta] = params.map(p => p.value);
			let S = 1; // Initial FX rate
			let sigma = Math.sqrt(omega / (1 - alpha - beta)); // Initial volatility (long-run mean)
			const dt = 0.01; // Time step
			const steps = 100; // Number of time steps
			const data = [S];

			for (let i = 0; i < steps; i++) {
				// Generate normal shock for return process
				const Z = gaussianRandom();

				// Update volatility process (conditional variance)
				const sigmaSquared = omega + alpha * (sigma * Z) ** 2 + beta * sigma ** 2;
				sigma = Math.sqrt(Math.max(sigmaSquared, 0)); // Ensure non-negative volatility

				// Update asset price
				const dS = sigma * Z * Math.sqrt(dt);
				S = S * Math.exp(dS); // Exponential model for price

				data.push(S);
			}

			return data;
		}
		// Simulate Variance Gamma Model
		function simulateVarianceGamma(params) {
			const [mu, sigma, theta, nu] = params.map(p => p.value);
			let S = 1; // Initial FX rate
			const dt = 0.01; // Time step
			const steps = 100; // Number of time steps
			const data = [S];

			for (let i = 0; i < steps; i++) {
				// Generate Gamma-distributed increment for time change
				const G_nu = gammaRandom(nu, dt);

				// Brownian increment
				const W_nu = Math.sqrt(dt) * gaussianRandom();

				// Update the asset price
				const driftTerm = (mu - 0.5 * sigma ** 2 - theta) * dt;
				const diffusionTerm = sigma * W_nu;
				const jumpTerm = theta * G_nu;

				S = S * Math.exp(driftTerm + diffusionTerm + jumpTerm);
				data.push(S);
			}

			return data;
		}


		// Gamma random number generator for the Gamma process
		function gammaRandom(nu, dt) {
			const shape = dt / nu;
			const scale = nu;
			let dGamma = 0;
			
			for (let i = 0; i < shape; i++) {
				dGamma += -Math.log(Math.random()) * scale;
			}
			
			return dGamma;
		}

        // Simulate Dupire Model (Local Volatility)
        function simulateDupire(params) {
            const [σ] = params.map(p => p.value);
            let S = 1; // Initial FX rate
            const dt = 0.01;
            const steps = 100;
            const data = [S];
            for (let i = 0; i < steps; i++) {
                const dW = Math.sqrt(dt) * gaussianRandom();
                S = S * Math.exp(-0.5 * σ ** 2 * dt + σ * dW);
                data.push(S);
            }
            return data;
        }
		
		function simulateMertonJump(params) {
			const [mu, sigma, lambda, nu, delta] = params.map(p => p.value);
			let S = 1; // Initial FX rate
			const dt = 0.01; // Time step
			const steps = 100; // Number of time steps
			const data = [S];

			for (let i = 0; i < steps; i++) {
				// Continuous Diffusion
				const dW = Math.sqrt(dt) * gaussianRandom();
				const diffusionPart = (mu - 0.5 * sigma ** 2) * dt + sigma * dW;

				// Jump Diffusion
				const dN = Math.random() < lambda * dt ? 1 : 0; // Poisson process (Bernoulli trial for jumps)
				const J = dN ? Math.exp(nu + delta * gaussianRandom()) : 1; // Jump size from log-normal if jump occurs

				// Update the asset price
				S = S * Math.exp(diffusionPart) * J;
				data.push(S);
			}

			return data;
		}

        // Simulate Extended CIR (ECIR) Model
        function simulateECIR(params) {
            const [κ, θ, σ] = params.map(p => p.value);
            let S = 1;
            let v = θ;
            const dt = 0.01;
            const steps = 100;
            const data = [S];
            for (let i = 0; i < steps; i++) {
                const dW = Math.sqrt(dt) * gaussianRandom();
                const dv = κ * (θ - v) * dt + σ * Math.sqrt(v) * dW;
                v = Math.max(v + dv, 0);
                S = S * Math.exp((v - 0.5 * v * dt) + Math.sqrt(v) * dW);
                data.push(S);
            }
            return data;
        }
				// Simulate Bates Model (Heston + Jumps)
		function simulateBates(params) {
			const [κ, θ, σ, ρ, v0, λ, ν, δ] = params.map(p => p.value);
			let S = 1;  // Initial FX rate
			let v = v0; // Initial variance
			const dt = 0.01; // Time step
			const steps = 100; // Number of time steps
			const data = [S];

			for (let i = 0; i < steps; i++) {
				const dW_S = Math.sqrt(dt) * gaussianRandom();
				const dW_v = Math.sqrt(dt) * gaussianRandom();
				const dW_corr = ρ * dW_S + Math.sqrt(1 - ρ**2) * dW_v; // Correlated Brownian motions

				// Variance (Heston)
				const dv = κ * (θ - v) * dt + σ * Math.sqrt(Math.max(v, 0)) * dW_corr;
				v = Math.max(v + dv, 0); // Variance cannot go below zero

				// Jump Diffusion
				const dN = Math.random() < λ * dt ? 1 : 0; // Poisson process (Bernoulli trial for jumps)
				const J =  Math.exp(ν + δ * gaussianRandom()); // Jump size from log-normal distribution

				// Asset price
				S = S * Math.exp((v - 0.5 * v * dt) + Math.sqrt(v) * dW_S) * (dN ? J : 1); // Apply jump if dN=1
				data.push(S);
			}

			return data;
		}


        // Gaussian random number generator for Brownian motion
        function gaussianRandom() {
            let u = Math.random(), v = Math.random();
            return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
        }

        // Initialize the app
        document.getElementById('model-select').addEventListener('change', updateModel);
        updateModel(); // Initial load
