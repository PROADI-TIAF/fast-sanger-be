/* Runs the Fast-Sanger-BE Python package inside a Web Worker with Pyodide.
 *
 * Messages in:  {type: "init"} | {type: "run", files: [{name, data}], guides, settings, names}
 * Messages out: {type: "status", text} | {type: "ready", version} | {type: "log", level, text}
 *               {type: "done", summary, zip, figures} | {type: "error", message}
 */
const PYODIDE_VERSION = "314.0.7";
const PYODIDE_URL = `https://cdn.jsdelivr.net/pyodide/v${PYODIDE_VERSION}/full/`;
const DISPLAYABLE = /\.(png|svg)$/i;

importScripts(PYODIDE_URL + "pyodide.js");

let booting = null;

function post(type, data = {}, transfer = []) {
  self.postMessage({ type, ...data }, transfer);
}

function pythonMessage(err) {
  // Keep the last line of a Python traceback ("WebInputError: ...") and drop the exception name.
  const lines = String(err && err.message ? err.message : err).trim().split("\n");
  return lines[lines.length - 1].replace(/^[\w.]*(Error|Exception): /, "");
}

async function boot() {
  post("status", { text: "Downloading Python (first visit only)" });
  const pyodide = await loadPyodide({ indexURL: PYODIDE_URL });
  post("status", { text: "Loading numpy, matplotlib and Biopython" });
  await pyodide.loadPackage(["numpy", "matplotlib", "biopython", "micropip"]);
  post("status", { text: "Installing Fast-Sanger-BE" });
  const manifest = await (await fetch("wheels/manifest.json", { cache: "no-cache" })).json();
  const urls = manifest.wheels.map((w) => new URL("wheels/" + w, self.location.href).href);
  await pyodide.runPythonAsync(
    `import json, micropip\nawait micropip.install(json.loads(${JSON.stringify(JSON.stringify(urls))}), deps=False)`
  );
  pyodide.runPython("import fast_sanger_be.web");
  const version = pyodide.runPython("import fast_sanger_be; fast_sanger_be.__version__");
  post("ready", { version });
  return pyodide;
}

function run(pyodide, msg) {
  pyodide.runPython("import shutil; shutil.rmtree('/work', ignore_errors=True)");
  pyodide.FS.mkdirTree("/work/input");
  for (const f of msg.files) pyodide.FS.writeFile(`/work/input/${f.name}`, new Uint8Array(f.data));
  pyodide.FS.writeFile("/work/input/guides.txt", msg.guides);

  const web = pyodide.pyimport("fast_sanger_be.web");
  let summary;
  try {
    const json = web.run(
      "/work/input", "/work/Final_Results", JSON.stringify(msg.settings), JSON.stringify(msg.names),
      (level, text) => post("log", { level, text })
    );
    summary = JSON.parse(json);
    web.zip_folder("/work/Final_Results", "/work/Final_Results.zip");
  } finally {
    web.destroy();
  }

  const zip = pyodide.FS.readFile("/work/Final_Results.zip");
  const figures = {};
  for (const blocks of Object.values(summary.editors)) {
    for (const block of blocks) {
      const paths = block.reads.map((r) => r.figure).concat(block.heatmap);
      for (const path of paths) {
        if (path && DISPLAYABLE.test(path)) figures[path] = pyodide.FS.readFile(`/work/Final_Results/${path}`).buffer;
      }
    }
  }
  post("done", { summary, zip: zip.buffer, figures }, [zip.buffer, ...Object.values(figures)]);
}

self.onmessage = async (event) => {
  const msg = event.data;
  let pyodide;
  try {
    if (!booting) booting = boot();
    pyodide = await booting;
  } catch (err) {
    booting = null; // allow a retry, e.g. after the connection comes back
    post("error", { message: "Could not load the analysis engine: " + pythonMessage(err), during: "boot" });
    return;
  }
  if (msg.type !== "run") return;
  try {
    run(pyodide, msg);
  } catch (err) {
    post("error", { message: pythonMessage(err), during: "run" });
  }
};
