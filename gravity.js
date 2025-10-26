const GravityDemo = {

  /*
  Basic simulation of a set of self-attracting particles which also repel as they get too close.
  The particles are drawn as discs with a changeable radius (see code).
  Two particles repel each other if their discs overlap.
  Two particles always attract one another regardless of distance between them.
  It is possible to have the attractive force decay as 1/R or 1/R/R (swappable interactively)
  */

  NAME: "Gravity",

  NPARTICLES: 30,
  RDISC: 0.400,

  GCONSTANT: 0.05,
  KCONSTANT: 20.0,
  MCONSTANT: 1.0,
  NUCONSTANT: 1.0e-3,

  LOGPOTENTIAL: true,
  USE_EULER: false,

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

    // Temporary storage required for predict-correct integrator
    this.VX_ = new Float64Array(this.NPARTICLES);
    this.VY_ = new Float64Array(this.NPARTICLES);
    this.FX_ = new Float64Array(this.NPARTICLES);
    this.FY_ = new Float64Array(this.NPARTICLES);

    this.reset_particle_state();
  },

  reset: function () {
    this.reset_particle_state();
  },

  evolve: function (t, dt) {
    if (this.USE_EULER) {
      this.calc_force(t, this.X, this.Y, this.VX, this.VY, this.FX, this.FY);
      for (var i = 0; i < this.NPARTICLES; i++) {
        this.X[i] += dt * this.VX[i];
        this.Y[i] += dt * this.VY[i];
        this.VX[i] += dt * this.FX[i] / this.MCONSTANT;
        this.VY[i] += dt * this.FY[i] / this.MCONSTANT;
      }
      return;
    }

    /*
    Predictor-corrector integrator examples

    * 2 evals per step (GPUSPH documentation):
    uh = u + (dt/2) * D(u)
    u  = u + dt * D(uh)
       = uh - (dt/2) * D(u) + dt * D(uh)
       = uh + dt * (D(uh) - 0.5 * D(u))

    * 2 evals per step (Heuns method)
    uf = u + dt * D(u)
    u  = u + (dt/2) * (D(u) + D(uf)) 
       = uf + (dt/2) * (D(uf) - D(u))
    */
    this.calc_force(t, this.X, this.Y, this.VX, this.VY, this.FX_, this.FY_);
    const half_dt = dt / 2.0;
    for (var i = 0; i < this.NPARTICLES; i++) {
      this.VX_[i] = this.VX[i];
      this.VY_[i] = this.VY[i];
      this.X[i] += half_dt * this.VX[i];
      this.Y[i] += half_dt * this.VY[i];
      this.VX[i] += half_dt * this.FX_[i] / this.MCONSTANT;
      this.VY[i] += half_dt * this.FY_[i] / this.MCONSTANT;
    }
    this.calc_force(t, this.X, this.Y, this.VX, this.VY, this.FX, this.FY);
    for (var i = 0; i < this.NPARTICLES; i++) {
      this.X[i] += dt * (this.VX[i] - 0.5 * this.VX_[i]);
      this.Y[i] += dt * (this.VY[i] - 0.5 * this.VY_[i]);
      this.VX[i] += dt * (this.FX[i] - 0.5 * this.FX_[i]) / this.MCONSTANT;
      this.VY[i] += dt * (this.FY[i] - 0.5 * this.FY_[i]) / this.MCONSTANT;
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
    // ctx.fillText("demo: " + this.NAME, 5, 35);
    ctx.fillText("[p] long-range force ~ " + (this.LOGPOTENTIAL ? "1/r" : "1/r/r") +
      ", [g|G] attract-strength: " + this.GCONSTANT.toFixed(4), 5, 35);
    ctx.fillText("[f|F] friction: " + this.NUCONSTANT.toFixed(6) +
      ", [k|K] repel-strength: " + this.KCONSTANT.toFixed(4), 5, 55);
    ctx.fillText("[e] integrator: " + (this.USE_EULER ? "Euler" : "predict-correct") +
      ", [d|D] disc radius: " + this.RDISC.toFixed(4), 5, 75);
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

    if (key == 'e' || key == 'E') {
      this.USE_EULER = !this.USE_EULER;
      return;
    }

    if (key == 'w' || key == 'W') {
      this.add_random_velocity();
      return;
    }

    if (key == '>') {
      this.add_angular_momentum(0.05);
      return;
    }

    if (key == '<') {
      this.add_angular_momentum(-0.05);
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

    if (key == 'g') {
      this.GCONSTANT *= 1.10;
      return;
    }

    if (key == 'G') {
      this.GCONSTANT /= 1.10;
      return;
    }

    if (key == 'k') {
      this.KCONSTANT *= 1.10;
      return;
    }

    if (key == 'K') {
      this.KCONSTANT /= 1.10;
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

  add_angular_momentum: function (omega) {
    var summ = 0.0;
    var sumx = 0.0;
    var sumy = 0.0;
    for (var i = 0; i < this.NPARTICLES; i++) {
      summ += this.MCONSTANT;
      sumx += this.X[i] * this.MCONSTANT;
      sumy += this.Y[i] * this.MCONSTANT;
    }
    const cgx = sumx / summ;
    const cgy = sumy / summ;
    for (var i = 0; i < this.NPARTICLES; i++) {
      const rx = this.X[i] - cgx;
      const ry = this.Y[i] - cgy;
      this.VX[i] -= omega * ry;
      this.VY[i] += omega * rx;
    }
  },

}
