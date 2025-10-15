/*
 * Basic demonstration of bouyancy calculations.
 *
 * Stokes wave representing a sea surface -- mocked-up perturbed hydrostatic pressure below surface. 
 * General quadrilateral cross-section of a floating object -- dynamic bouyancy (boundary integral).
 *
 */

// TODO: this set of codes should be packaged in a namespace-like object & expose a set of callbacks
//       which include, init, draw, evove, draw_stats, handle_keypress, and such...

/*
const MyNamespace = {
  myFunction: () => console.log("Hello!"),
  myVariable: 42
};

MyNamespace.myFunction(); // Output: Hello!
console.log(MyNamespace.myVariable); // Output: 42
*/

const PRES_SURFACE = 101.0e3; // arbitrary
const GRAVITY_ACC = 9.82;
const ATMOSPHERE_RHO = 1.225;
const SEAWATER_RHO = 1025.0;
const OAK_RHO = 700.0;
const REDWOOD_RHO = 450.0;
const BALSA_RHO = 150.0;

let WAVE_LAMBDA = 10.0;
let WAVE_K = 2 * Math.PI / WAVE_LAMBDA;
let WAVE_A = (0.5 * WAVE_LAMBDA * 0.142) / 5; // replace the 3 with 1 to get "maxed out wave amplitude"
let WAVE_PTS = 257;

function stokes_3rd(ka, theta) {
    const one_ = Math.cos(theta) * (1.0 - (ka * ka) / 16.0);
    const two_ = Math.cos(2 * theta) * 0.5 * ka;
    const three_ = Math.cos(3 * theta) * (3.0 / 8.0) * ka * ka;
    return one_ + two_ + three_;
}

function stokes_omega(k, a) {
    const ka = k * a;
    const c = (1.0 + 0.5 * ka * ka) * Math.sqrt(GRAVITY_ACC / k); // c = omega / k
    return c / k;
}

function perturbed_pressure(x, y, k, a, omega, t) {
    // Express the pressure in the (x,y) plane given the pertubed sea surface by the Stokes wave.
    // Not an exact expression but good enough for a nice looking animation.
    // Ignore the atmospheric pressure, let it be zero for simplicity.
    const theta = k * x - omega * t;
    const eta = a * stokes_3rd(k * a, theta);
    if (y >= eta) {
        const h = y - eta;
        const eff_h = y - eta * Math.exp(-k * h);
        return -1 * eff_h * ATMOSPHERE_RHO * GRAVITY_ACC + PRES_SURFACE;
    }
    const h = eta - y;
    const eff_h = eta * Math.exp(-k * h) - y;
    return eff_h * SEAWATER_RHO * GRAVITY_ACC + PRES_SURFACE;
}

function draw_stokes_wave(t, k, a, ctx, xmin, xmax, npts) {
    ctx.strokeStyle = "blue";
    ctx.lineWidth = 0.020;
    const omega = stokes_omega(k, a);
    const ka = k * a;
    const dx = (xmax - xmin) / (npts - 1);
    ctx.beginPath();
    ctx.moveTo(xmin, a * stokes_3rd(ka, k * xmin - omega * t));
    for (var i = 1; i < npts; i++) {
        const xi = xmin + dx * i;
        ctx.lineTo(xi, a * stokes_3rd(ka, k * xi - omega * t));
    }
    ctx.stroke();
    ctx.closePath();
}

//
// TODO: Clenshaw curtis quadrature for better accuracy (?!)
//

// Basic quadrature nodes/weights on interval [-1, 1]; n >= 2
function trapezoidal_nodes_and_weights(n) {
    const x = Array(n).fill(0.0);
    const dx = 2.0 / (x.length - 1);
    for (var i = 0; i < x.length; i++) {
        x[i] = -1.0 + dx * i;
    }
    const w = Array(n).fill(dx);
    w[0] = 0.5 * dx;
    w[w.length - 1] = 0.5 * dx;
    return [x, w];
}

function single_line_integral(ax, ay, bx, by, nodes, weights, pres_params) {
    // pres_params = [k, a, omega, t]
    var A = 0.0;
    var AX = 0.0;
    var AY = 0.0;
    var AXX = 0.0;
    var AYY = 0.0;
    var L = 0.0;
    var P = 0.0;
    var PX = 0.0;
    var PY = 0.0;
    var WL = 0.0;
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
        const p_ = perturbed_pressure(x_, y_, pres_params[0], pres_params[1], pres_params[2], pres_params[3]);
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
        WL += wdl_ * (p_ > PRES_SURFACE ? 1.0 : 0.0); // wet perimeter
    }
    return [A, AX, AY, AXX, AYY, L, P * nxhat, P * nyhat, PX * nyhat, PY * nxhat, WL];
}

// Line-integration around the boundary af the quadrilateral ABCDA
// (x1,w1) are quadrature nodes/weights (reused for all 4 boundary lines)
function ql_quadrature(qlx, qly, x1, w1, pres_params) {
    if (qlx.length != 4 || qlx.length != qly.length) return [];
    if (x1.length != w1.length) return [];

    // A->B
    const AB = single_line_integral(qlx[0], qly[0], qlx[1], qly[1], x1, w1, pres_params);

    // B->C
    const BC = single_line_integral(qlx[1], qly[1], qlx[2], qly[2], x1, w1, pres_params);

    // C->D
    const CD = single_line_integral(qlx[2], qly[2], qlx[3], qly[3], x1, w1, pres_params);

    // D->A
    const DA = single_line_integral(qlx[3], qly[3], qlx[0], qly[0], x1, w1, pres_params);

    var Q = Array(AB.length).fill(0.0);
    for (var i = 0; i < AB.length; i++) {
        Q[i] = AB[i] + BC[i] + CD[i] + DA[i];
    }

    return Q;
}


function createQuadrilateralObject() {
    let obj = {
        // Cross-section of object is a general quadrilateral qlx/y = x/y coordinates
        // Specify vertices in CCW order to get the boundary integration sign correct
        qlx: [-1.0, 1.0, 1.0, -1.0],
        qly: [-0.5, -0.5, 0.5, 0.5],

        // uniform density of material
        rho: REDWOOD_RHO,

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

        quadrature_arrays: trapezoidal_nodes_and_weights(100),

        update_mechanics: function (time, pres_params) {
            const line_integrals = ql_quadrature(this.qlx,
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
            for (var i = 0; i < 4; i++) {
                const dx_ = this.qlx[i] - this.cgx;
                const dy_ = this.qly[i] - this.cgy;
                const x_ = costh * dx_ + sinth * dy_;
                const y_ = -sinth * dx_ + costh * dy_;
                this.qlx[i] = x_ + this.cgx + delta_time * this.vx;
                this.qly[i] = y_ + this.cgy + delta_time * this.vy;
            }

            const nu_frac = this.wet_perimeter / this.perimeter;

            this.vx += delta_time * (this.pres_Fx - this.vx * this.nu * nu_frac) / this.mass;
            this.vy += delta_time * (this.pres_Fy - this.mass * GRAVITY_ACC - this.vy * this.nu * nu_frac) / this.mass;
            this.omegaz += delta_time * (this.pres_Tz - this.omegaz * this.nuz * nu_frac) / this.Iz;
        },
    };
    return obj;
}

function drawQuadrilateralObject(obj, ctx) {
    ctx.strokeStyle = "rgba(255,255,255,0.80)";
    ctx.lineWidth = 0.040;
    ctx.beginPath();
    ctx.moveTo(obj.qlx[0], obj.qly[0]);
    ctx.lineTo(obj.qlx[1], obj.qly[1]);
    ctx.lineTo(obj.qlx[2], obj.qly[2]);
    ctx.lineTo(obj.qlx[3], obj.qly[3]);
    ctx.lineTo(obj.qlx[0], obj.qly[0]);
    ctx.stroke();
    ctx.closePath();

    ctx.fillStyle = "rgba(255,255,255,0.50)";
    ctx.beginPath();
    ctx.arc(obj.cgx, obj.cgy, 0.10, 0, 2 * Math.PI);
    ctx.closePath();
    ctx.fill();

    ctx.fillStyle = "rgba(128,128,255,0.50)";
    ctx.beginPath();
    ctx.arc(obj.cbx, obj.cby, 0.10, 0, 2 * Math.PI);
    ctx.closePath();
    ctx.fill();
}
