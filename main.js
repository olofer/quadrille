/*

Demo of the planeview.js grid plotter with an interactive animation.

...

TODO: ability to change the side length of the quadrilateral
TODO: ability to change wave amp and wavelength (with safeguards)

*/

var floater = createQuadrilateralObject();
floater.update_mechanics(0.0, [WAVE_K, WAVE_A, stokes_omega(WAVE_K, WAVE_A), 0.0]);

console.log(floater);
console.log(floater.rho * 5 * floater.area / 12); // NOTE: should (ideally) match floater.Iz (for default quadrilateral)

function printSimulatorStats(ctx, timestamp, fpsval) {
    ctx.font = "17px Arial";
    ctx.fillStyle = "white";
    ctx.fillText("time: " + timestamp.toFixed(3) + " sec" + " [" + fpsval.toFixed(1) + " fps]", 5, 15);

    ctx.fillText("mass: " + floater.mass.toFixed(3) +
        " [kg], cg = (" + floater.cgx.toFixed(3) +
        "," + floater.cgy.toFixed(3) + ")", 5, 35);

    ctx.fillText("Iz[cg]: " + floater.Iz.toFixed(3) + ", area: " + floater.area.toFixed(3) + " [m^2]", 5, 55);

    ctx.fillText("perimeter (of which wet): " + floater.perimeter.toFixed(3) +
        " (" + floater.wet_perimeter.toFixed(3) + ") [m]", 5, 75);
}

const DEFAULT_BBOX = [-5.0, 5.0, -3.0, 3.0];

const canvas = document.getElementById("canvas");
var PV = createDefaultPlaneView(canvas.width, canvas.height);

var viewZoom = 1.00;
var viewEta = 0.125;
var bbox = [DEFAULT_BBOX[0], DEFAULT_BBOX[1], DEFAULT_BBOX[2], DEFAULT_BBOX[3]];

let FPS = 50.0;  // 20ms refresh intervals
let dt = 0.0020;  // 2.0ms simulator time-stepping

let filtered_fpsval = 0.0;
const filter_beta = 0.990;

let tsim = 0.0;
let frame = 0;
var startTime = Date.now();

function resizeWindow() {
    canvas.width = window.innerWidth - 24;
    canvas.height = window.innerHeight - 24;
    PV = createDefaultPlaneView(canvas.width, canvas.height);
}

function reset_state() {
    tsim = 0.0;
    frame = 0;
    startTime = Date.now();
}

function refresh() {

    var currentTime = Date.now();
    var elapsedTime = currentTime - startTime;
    startTime = currentTime;

    if (elapsedTime == 0.0) return;

    var elapsedSec = elapsedTime / 1000.0;

    while (elapsedSec > 0.0) {

        floater.update_mechanics(tsim, [WAVE_K, WAVE_A, stokes_omega(WAVE_K, WAVE_A), tsim]);
        floater.evolve(tsim, dt);

        tsim += dt;
        elapsedSec -= dt;
    }

    PV.autoZoom(bbox[0], bbox[1], bbox[2], bbox[3], viewZoom, viewEta);

    const ctx = canvas.getContext("2d");

    PV.drawGrid(ctx);
    PV.setTransform(ctx);

    draw_stokes_wave(tsim, WAVE_K, WAVE_A, ctx, PV.xmin, PV.xmax, WAVE_PTS);
    drawQuadrilateralObject(floater, ctx);

    PV.unitTransform(ctx);
    filtered_fpsval = filter_beta * filtered_fpsval + (1.0 - filter_beta) * (1000.0 / elapsedTime);
    printSimulatorStats(ctx, tsim, filtered_fpsval);
    PV.setTransform(ctx);

    frame++;
}

function keyDownEvent(e) {
    var code = e.keyCode;
    var key = e.key;

    if (key == 'r' || key == 'R') {
        reset_state();
        return;
    }

    if (key == 'z' || key == 'Z') { // reset view
        bbox = [DEFAULT_BBOX[0], DEFAULT_BBOX[1], DEFAULT_BBOX[2], DEFAULT_BBOX[3]];
        viewZoom = 1.0;
        return;
    }

    if (key == 'c' || key == 'C') { // center view (do not reset zoom level)
        bbox = [DEFAULT_BBOX[0], DEFAULT_BBOX[1], DEFAULT_BBOX[2], DEFAULT_BBOX[3]];
        return;
    }

    if (key == 'm' || key == 'M') { // swap type of wood
        floater.rho = floater.rho == WOOD_RHO ? BALSA_RHO : WOOD_RHO;
        return;
    }

    if (key == 'a') {
        WAVE_A *= 0.90;
        return;
    }

    if (key == 'A') {
        WAVE_A *= 1.10;
        return;
    }

    if (key == ' ' && !e.shiftKey) {
        floater.omegaz += 2.0;
        return;
    }

    if (key == ' ' && e.shiftKey) {
        floater.omegaz -= 2.0;
        return;
    }

    if (key == '>') {
        floater.vx += 1.0;
        return;
    }

    if (key == '<') {
        floater.vx -= 1.0;
        return;
    }

    if (code == 38 && e.shiftKey) {
        viewZoom *= 0.80;
        return;
    }

    if (code == 40 && e.shiftKey) {
        viewZoom *= 1.25;
        return;
    }

    const dx = (bbox[1] - bbox[0]) / 12.0;
    const dy = (bbox[3] - bbox[2]) / 12.0;
    const dmove = viewZoom * (dx + dy) / 2;

    if (code === 39) // right
    {
        bbox[0] += dmove;
        bbox[1] += dmove;
    }
    else if (code === 37) // left
    {
        bbox[0] -= dmove;
        bbox[1] -= dmove;
    }
    else if (code === 38) // up
    {
        bbox[2] += dmove;
        bbox[3] += dmove;
    }
    else if (code === 40) // down
    {
        bbox[2] -= dmove;
        bbox[3] -= dmove;
    }
}

function keyUpEvent(e) {
    var code = e.keyCode;
    var key = e.key;
}

resizeWindow();
reset_state();

window.addEventListener('resize', resizeWindow);
window.addEventListener('keydown', keyDownEvent);
window.addEventListener('keyup', keyUpEvent);

let refresher_id = setInterval(refresh, 1000 / FPS);
