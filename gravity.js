const GravityDemo = {

  /*
  This is a basic simulation of a set of self-attracting particles which also repel as they get too close.

  TODO: key press "<" and ">" to impart angular momentum around the present center of gravity of the system
  TODO: more reliable integrator (predictor-corrector type or just explicit Heun)
  TODO: ability to change  G strength and K strength
  */

  NAME: "Gravity",

  NPARTICLES: 30,
  RDISC: 0.400,

  GCONSTANT: 0.05,
  KCONSTANT: 20.0,
  MCONSTANT: 1.0,
  NUCONSTANT: 1.0e-3,

  LOGPOTENTIAL: true,

  X: null,
  Y: null,
  VX: null,
  VY: null,
  FX: null,
  FY: null,

  init: function () {
    console.log("Initializing: " + this.NAME);

    this.X = new Float64Array(this.NPARTICLES);
    this.Y = new Float64Array(this.NPARTICLES);
    this.VX = new Float64Array(this.NPARTICLES);
    this.VY = new Float64Array(this.NPARTICLES);
    this.FX = new Float64Array(this.NPARTICLES);
    this.FY = new Float64Array(this.NPARTICLES);

    this.reset_particle_state();
  },

  reset: function () {
    this.reset_particle_state();
  },

  evolve: function (t, dt) {
    this.calc_force(t, this.X, this.Y, this.VX, this.VY, this.FX, this.FY);
    for (var i = 0; i < this.NPARTICLES; i++) {
      this.X[i] += dt * this.VX[i];
      this.Y[i] += dt * this.VY[i];
      this.VX[i] += dt * this.FX[i] / this.MCONSTANT;
      this.VY[i] += dt * this.FY[i] / this.MCONSTANT;
    }
  },

  draw: function (ctx, t, xmin, xmax, ymin, ymax) {

    ctx.fillStyle = `rgba(200, 255, 200, 0.67)`; // red with opacity
    for (var i = 0; i < this.NPARTICLES; i++) {
      ctx.beginPath();
      ctx.arc(this.X[i], this.Y[i], this.RDISC, 0, 2 * Math.PI);
      ctx.fill();
    }
  },

  print_stats: function (ctx, t) {
    ctx.fillText("demo: " + this.NAME, 5, 35);
    ctx.fillText("[p] long-range force ~ " + (this.LOGPOTENTIAL ? "1/r" : "1/r/r"), 5, 55);
    ctx.fillText("[f] fricton: " + this.NUCONSTANT.toFixed(6), 5, 75);
  },

  handle_key_down: function (e) {
    const code = e.keyCode;
    const key = e.key;

    if (key == ' ' && !e.shiftKey) {
      this.zero_velocity();
      return;
    }

    if (key == 'p' || key == 'P') {
      this.LOGPOTENTIAL = !this.LOGPOTENTIAL;
      return;
    }

    if (key == 'w' || key == 'W') {
      this.add_random_velocity();
      return;
    }

    if (key == 'd') {
      this.RDISC *= 1.10;
      return;
    }

    if (key == 'D') {
      this.RDISC /= 1.10;
      return;
    }

    if (key == 'f') {
      this.NUCONSTANT *= 2.0;
      return;
    }

    if (key == 'F') {
      this.NUCONSTANT /= 2.0;
      return;
    }
  },

  calc_force: function (t, X, Y, VX, VY, FX, FY) {
    const CONSTANT = this.GCONSTANT * this.MCONSTANT * this.MCONSTANT;
    const EPS = 0.001;
    const KAPPA = (EPS + this.RDISC) / (EPS + this.RDISC * this.RDISC);
    var Fij = 0.0;
    for (var i = 0; i < this.NPARTICLES; i++) {
      FX[i] = 0.0;
      FY[i] = 0.0;
      const xi = X[i];
      const yi = Y[i];
      for (var j = 0; j < this.NPARTICLES; j++) {
        if (j == i) continue;
        const xj = X[j];
        const yj = Y[j];
        const dx = xj - xi;
        const dy = yj - yi;
        const dsq = dx * dx + dy * dy;
        const d = Math.sqrt(dsq);
        const rhatx = dx / d;
        const rhaty = dy / d;
        if (!this.LOGPOTENTIAL)
          Fij = CONSTANT / (EPS + dsq);
        else
          Fij = CONSTANT * KAPPA / (EPS + d);
        FX[i] += Fij * rhatx;
        FY[i] += Fij * rhaty;
        if (d < 2 * this.RDISC) {
          const z = 2 * this.RDISC - d;
          FX[i] -= this.KCONSTANT * z * rhatx;
          FY[i] -= this.KCONSTANT * z * rhaty;
        }
        FX[i] -= this.NUCONSTANT * VX[i];
        FY[i] -= this.NUCONSTANT * VY[i];
      }
    }
  },

  reset_particle_state: function () {
    // TODO: this one need to set up VX, VY so that there is an angular momentum to keep things more tidy..
    for (var i = 0; i < this.NPARTICLES; i++) {
      this.X[i] = (2.0 * Math.random() - 1.0) * 5.0;
      this.Y[i] = (2.0 * Math.random() - 1.0) * 5.0;
      this.VX[i] = (2.0 * Math.random() - 1.0) * 0.5;
      this.VY[i] = (2.0 * Math.random() - 1.0) * 0.5;
      this.FX[i] = 0.0;
      this.FY[i] = 0.0;
    }
  },

  zero_velocity: function () {
    for (var i = 0; i < this.NPARTICLES; i++) {
      this.VX[i] = 0.0;
      this.VY[i] = 0.0;
      this.FX[i] = 0.0;
      this.FY[i] = 0.0;
    }
  },

  add_random_velocity: function () {
    for (var i = 0; i < this.NPARTICLES; i++) {
      this.VX[i] += (2.0 * Math.random() - 1.0) * 0.5;
      this.VY[i] += (2.0 * Math.random() - 1.0) * 0.5;
      this.FX[i] = 0.0;
      this.FY[i] = 0.0;
    }
  },

}
