/*
 * Basic demonstration of bouyancy calculations.
 *
 * Stokes wave representing a sea surface -- mocked-up perturbed hydrostatic pressure below surface. 
 * General quadrilateral cross-section of a floating object -- dynamic bouyancy (boundary integral).
 *
 */

function createQuadrilateralObject(demo_object) {
    let obj = {
        // Cross-section of object is a general quadrilateral qlx/y = x/y coordinates
        // Specify vertices in CCW order to get the boundary integration sign correct
        qlx: [-1.0, 1.0, 1.0, -1.0],
        qly: [-0.5, -0.5, 0.5, 0.5],

        // uniform density of material
        rho: demo_object.REDWOOD_RHO,

        area: 0.0, // cross section area
        perimeter: 0.0,
        mass: 0.0, // mass (per meter into plane)
        wet_perimeter: 0.0,

        // center of gravity
        cgx: 0.0,
        cgy: 0.0,

        pres_Fx: 0.0,
        pres_Fy: 0.0,
        pres_Tz: 0.0,

        // center of bouyancy
        cbx: 0.0,
        cby: 0.0,

        // moment of inertia (around CG)
        Iz: 0.0,

        // damping parameter (full utilization if the object is fully underwater)
        nu: 5000.00,
        nuz: 5000.00,

        // angular and linear velocities
        omegaz: 1.55,
        vx: -0.10,
        vy: 5.00,

        quadrature_arrays: demo_object.trapezoidal_nodes_and_weights(100),

        update_mechanics: function (time, pres_params) {
            const line_integrals = demo_object.ql_quadrature(this.qlx,
                this.qly,
                this.quadrature_arrays[0],
                this.quadrature_arrays[1],
                pres_params);

            // line_integrals = [A, AX, AY, AXX, AYY, L, ...]

            this.area = line_integrals[0];
            this.mass = this.area * this.rho;
            this.cgx = line_integrals[1] / this.area;
            this.cgy = line_integrals[2] / this.area;

            const Azz_ = line_integrals[3] + line_integrals[4]; // polar area moment
            const dsq_ = this.cgx * this.cgx + this.cgy * this.cgy;

            // Parallel axis translation to get CG moment of inertia
            this.Iz = (Azz_ - dsq_ * this.area) * this.rho;

            this.perimeter = line_integrals[5];
            this.wet_perimeter = line_integrals[10];

            this.pres_Fx = line_integrals[6];
            this.pres_Fy = line_integrals[7];

            // Translate integrals to torque w.r.t. CM
            const part_1 = line_integrals[8] - this.cgx * line_integrals[7];
            const part_2 = line_integrals[9] - this.cgy * line_integrals[6];

            this.pres_Tz = -part_1 + part_2;

            // So that the CB can be visualized
            this.cbx = line_integrals[8] / this.pres_Fy;
            this.cby = line_integrals[9] / this.pres_Fx;
        },

        // Time-stepper not appropriate for accuracy, but fine for simple animations.
        evolve: function (time, delta_time) {
            const dtheta = this.omegaz * delta_time;
            const sinth = Math.sin(dtheta);
            const costh = Math.cos(dtheta);
            for (let i = 0; i < 4; i++) {
                const dx_ = this.qlx[i] - this.cgx;
                const dy_ = this.qly[i] - this.cgy;
                const x_ = costh * dx_ + sinth * dy_;
                const y_ = -sinth * dx_ + costh * dy_;
                this.qlx[i] = x_ + this.cgx + delta_time * this.vx;
                this.qly[i] = y_ + this.cgy + delta_time * this.vy;
            }

            const nu_frac = this.wet_perimeter / this.perimeter;

            this.vx += delta_time * (this.pres_Fx - this.vx * this.nu * nu_frac) / this.mass;
            this.vy += delta_time * (this.pres_Fy - this.mass * demo_object.GRAVITY_ACC - this.vy * this.nu * nu_frac) / this.mass;
            this.omegaz += delta_time * (this.pres_Tz - this.omegaz * this.nuz * nu_frac) / this.Iz;
        },

        draw: function (ctx) {
            ctx.strokeStyle = "rgba(255,255,255,0.80)";
            ctx.lineWidth = 0.040;
            ctx.beginPath();
            ctx.moveTo(this.qlx[0], this.qly[0]);
            ctx.lineTo(this.qlx[1], this.qly[1]);
            ctx.lineTo(this.qlx[2], this.qly[2]);
            ctx.lineTo(this.qlx[3], this.qly[3]);
            ctx.lineTo(this.qlx[0], this.qly[0]);
            ctx.stroke();
            ctx.closePath();

            ctx.fillStyle = "rgba(255,255,255,0.50)";
            ctx.beginPath();
            ctx.arc(this.cgx, this.cgy, 0.10, 0, 2 * Math.PI);
            ctx.closePath();
            ctx.fill();

            ctx.fillStyle = "rgba(128,128,255,0.50)";
            ctx.beginPath();
            ctx.arc(this.cbx, this.cby, 0.10, 0, 2 * Math.PI);
            ctx.closePath();
            ctx.fill();
        },
    };
    return obj;
}

// TODO: ability to change the side length of the quadrilateral
// TODO: ability to change wave amp and wavelength (with safeguards)

const BouyancyDemo = {

    NAME: "Bouyancy",

    PRES_SURFACE: 101.0e3, // arbitrary
    GRAVITY_ACC: 9.82,
    ATMOSPHERE_RHO: 1.225,
    SEAWATER_RHO: 1025.0,
    OAK_RHO: 700.0,
    REDWOOD_RHO: 450.0,
    BALSA_RHO: 150.0,

    WAVE_LAMBDA: 0.0,
    WAVE_K: 0.0,
    WAVE_A: 0.0,
    WAVE_PTS: 0,

    floater: null,

    init: function () {
        console.log("Initializing: " + this.NAME);

        this.WAVE_LAMBDA = 10.0;
        this.WAVE_K = 2 * Math.PI / this.WAVE_LAMBDA;
        this.WAVE_A = (0.5 * this.WAVE_LAMBDA * 0.142) / 5; // replace the 3 with 1 to get "maxed out wave amplitude"
        this.WAVE_PTS = 257;

        this.floater = createQuadrilateralObject(this);
        this.floater.update_mechanics(0.0, [this.WAVE_K, this.WAVE_A, this.stokes_omega(this.WAVE_K, this.WAVE_A), 0.0]);

        console.log(this.floater);

        // NOTE: should (ideally) match floater.Iz (for default quadrilateral)
        console.log(this.floater.rho * 5 * this.floater.area / 12);
    },

    reset: function () {
        // TODO: something?
    },

    evolve: function (tsim_, dt_) {
        this.floater.update_mechanics(tsim_, [this.WAVE_K, this.WAVE_A, this.stokes_omega(this.WAVE_K, this.WAVE_A), tsim_]);
        this.floater.evolve(tsim_, dt_);
    },

    draw: function (ctx_, tsim_, xmin_, xmax_, ymin_, ymax_) {
        this.draw_stokes_wave(tsim_, this.WAVE_K, this.WAVE_A, ctx_, xmin_, xmax_, this.WAVE_PTS);
        this.floater.draw(ctx_);
    },

    print_stats: function (ctx_, tsim_) {
        ctx_.fillText(
            "mass: " + this.floater.mass.toFixed(3) +
            " [kg], cg = (" + this.floater.cgx.toFixed(3) +
            "," + this.floater.cgy.toFixed(3) + ")", 5, 35);

        ctx_.fillText(
            "Iz[cg]: " + this.floater.Iz.toFixed(3) +
            ", area: " + this.floater.area.toFixed(3) + " [m^2]", 5, 55);

        ctx_.fillText(
            "perimeter (of which wet): " + this.floater.perimeter.toFixed(3) +
            " (" + this.floater.wet_perimeter.toFixed(3) + ") [m]", 5, 75);
    },

    handle_key_down: function (e) {
        const code = e.keyCode;
        const key = e.key;

        if (key == 'm' || key == 'M') { // swap type of wood
            if (this.floater.rho == this.BALSA_RHO)
                this.floater.rho = this.REDWOOD_RHO;
            else if (this.floater.rho == this.REDWOOD_RHO)
                this.floater.rho = this.OAK_RHO;
            else if (this.floater.rho == this.OAK_RHO)
                this.floater.rho = this.BALSA_RHO;
            return;
        }

        if (key == 'a') {
            this.WAVE_A *= 0.90;
            return;
        }

        if (key == 'A') {
            this.WAVE_A *= 1.10;
            return;
        }

        if (key == ' ' && !e.shiftKey) {
            this.floater.omegaz += 2.0;
            return;
        }

        if (key == ' ' && e.shiftKey) {
            this.floater.omegaz -= 2.0;
            return;
        }

        if (key == '>') {
            this.floater.vx += 1.0;
            return;
        }

        if (key == '<') {
            this.floater.vx -= 1.0;
            return;
        }
    },

    stokes_3rd: function (ka, theta) {
        const one_ = Math.cos(theta) * (1.0 - (ka * ka) / 16.0);
        const two_ = Math.cos(2 * theta) * 0.5 * ka;
        const three_ = Math.cos(3 * theta) * (3.0 / 8.0) * ka * ka;
        return one_ + two_ + three_;
    },

    stokes_omega: function (k, a) {
        const ka = k * a;
        const c = (1.0 + 0.5 * ka * ka) * Math.sqrt(this.GRAVITY_ACC / k); // c = omega / k
        return c / k;
    },

    perturbed_pressure: function (x, y, k, a, omega, t) {
        // Express the pressure in the (x,y) plane given the pertubed sea surface by the Stokes wave.
        // Not an exact expression but good enough for a nice looking animation.
        // Ignore the atmospheric pressure, let it be zero for simplicity.
        const theta = k * x - omega * t;
        const eta = a * this.stokes_3rd(k * a, theta);
        if (y >= eta) {
            const h = y - eta;
            const eff_h = y - eta * Math.exp(-k * h);
            return -1 * eff_h * this.ATMOSPHERE_RHO * this.GRAVITY_ACC + this.PRES_SURFACE;
        }
        const h = eta - y;
        const eff_h = eta * Math.exp(-k * h) - y;
        return eff_h * this.SEAWATER_RHO * this.GRAVITY_ACC + this.PRES_SURFACE;
    },

    draw_stokes_wave: function (t, k, a, ctx, xmin, xmax, npts) {
        ctx.strokeStyle = "blue";
        ctx.lineWidth = 0.020;
        const omega = this.stokes_omega(k, a);
        const ka = k * a;
        const dx = (xmax - xmin) / (npts - 1);
        ctx.beginPath();
        ctx.moveTo(xmin, a * this.stokes_3rd(ka, k * xmin - omega * t));
        for (let i = 1; i < npts; i++) {
            const xi = xmin + dx * i;
            ctx.lineTo(xi, a * this.stokes_3rd(ka, k * xi - omega * t));
        }
        ctx.stroke();
        ctx.closePath();
    },

    trapezoidal_nodes_and_weights: function (n) {
        // Basic quadrature nodes/weights on interval [-1, 1]; n >= 2
        const x = Array(n).fill(0.0);
        const dx = 2.0 / (x.length - 1);
        for (let i = 0; i < x.length; i++) {
            x[i] = -1.0 + dx * i;
        }
        const w = Array(n).fill(dx);
        w[0] = 0.5 * dx;
        w[w.length - 1] = 0.5 * dx;
        return [x, w];
    },

    single_line_integral: function (ax, ay, bx, by, nodes, weights, pres_params) {
        // pres_params = [k, a, omega, t]
        let A = 0.0;
        let AX = 0.0;
        let AY = 0.0;
        let AXX = 0.0;
        let AYY = 0.0;
        let L = 0.0;
        let P = 0.0;
        let PX = 0.0;
        let PY = 0.0;
        let WL = 0.0;
        const midx = (ax + bx) / 2;
        const midy = (ay + by) / 2;
        const tx = (bx - ax) / 2;
        const ty = (by - ay) / 2;
        const nx = -ty;
        const ny = tx;
        const dl = Math.sqrt(tx * tx + ty * ty); // length per unit of line parameter
        const nxhat = nx / dl;
        const nyhat = ny / dl;
        for (var i = 0; i < nodes.length; i++) {
            const t_ = nodes[i]; // line parameter t = -1 .. +1
            const x_ = midx + t_ * tx;
            const y_ = midy + t_ * ty;
            const p_ = this.perturbed_pressure(x_, y_, pres_params[0], pres_params[1], pres_params[2], pres_params[3]);
            const w_ = weights[i];
            const xy_ = x_ * y_;
            const xx_ = x_ * x_;
            const yy_ = y_ * y_;
            A += w_ * (ty * x_ - y_ * tx) / 2; // area integral piece
            AX += w_ * (-1 * xy_ * tx + 0.5 * xx_ * ty) / 2;
            AY += w_ * (xy_ * ty - tx * 0.5 * yy_) / 2;
            AXX += w_ * (-1 * xy_ * x_ * tx + ty * xx_ * x_ / 3) / 2;
            AYY += w_ * (xy_ * y_ * ty - tx * yy_ * y_ / 3) / 2;
            const wdl_ = w_ * dl;
            const wdlp_ = wdl_ * p_;
            L += wdl_;
            P += wdlp_;
            PX += wdlp_ * x_;
            PY += wdlp_ * y_;
            WL += wdl_ * (p_ > this.PRES_SURFACE ? 1.0 : 0.0); // wet perimeter
        }
        return [A, AX, AY, AXX, AYY, L, P * nxhat, P * nyhat, PX * nyhat, PY * nxhat, WL];
    },

    // Line-integration around the boundary af the quadrilateral ABCDA
    // (x1,w1) are quadrature nodes/weights (reused for all 4 boundary lines)
    ql_quadrature: function (qlx, qly, x1, w1, pres_params) {
        if (qlx.length != 4 || qlx.length != qly.length) return [];
        if (x1.length != w1.length) return [];

        // A->B
        const AB = this.single_line_integral(qlx[0], qly[0], qlx[1], qly[1], x1, w1, pres_params);

        // B->C
        const BC = this.single_line_integral(qlx[1], qly[1], qlx[2], qly[2], x1, w1, pres_params);

        // C->D
        const CD = this.single_line_integral(qlx[2], qly[2], qlx[3], qly[3], x1, w1, pres_params);

        // D->A
        const DA = this.single_line_integral(qlx[3], qly[3], qlx[0], qly[0], x1, w1, pres_params);

        var Q = Array(AB.length).fill(0.0);
        for (let i = 0; i < AB.length; i++) {
            Q[i] = AB[i] + BC[i] + CD[i] + DA[i];
        }

        return Q;
    },

}
