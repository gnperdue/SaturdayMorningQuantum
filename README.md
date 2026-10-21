# Quantum Computing in 2 Hours

Educational Jupyter notebooks teaching quantum computing concepts in a
workshop format. Topics covered include quantum mechanics foundations,
classical computing background, quantum gates, and fault-tolerant
quantum computing. There are two notebooks:

- **`QuantumComputing.ipynb`** — the lecture: OpenQASM programming on a
  simulator, digital (fault-tolerant) quantum computing and its
  algorithms, two demos (a quantum adder and a rotation built from H, S
  and T gates, with a precessing spin), and a bonus "Quantum Spy Hunter"
  key distribution application.
- **`QuantumComputing_QEC.ipynb`** — the companion on quantum error
  correction: repetition codes, the surface code and its threshold, and
  how magic states and T factories make fault-tolerant gates (uses Cirq,
  Stim and PyMatching).

## Prerequisites

- [uv](https://docs.astral.sh/uv/) (Python package and project manager)

That's it — `uv` handles Python and all dependencies automatically.

## Quick start

### 1. Install uv

macOS / Linux:

```sh
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Windows:

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

See the [uv installation docs](https://docs.astral.sh/uv/getting-started/installation/)
for other methods (Homebrew, pipx, etc.).

### 2. Clone the repository

```sh
git clone https://github.com/gnperdue/SaturdayMorningQuantum.git
cd SaturdayMorningQuantum
```

### 3. Launch JupyterLab

```sh
uv run jupyter lab QuantumComputing.ipynb
```

On the first run, `uv` will automatically:

1. Download Python 3.14 (if not already installed)
2. Create a virtual environment in `.venv/`
3. Install all dependencies from `pyproject.toml`
4. Start JupyterLab

Open `QuantumComputing.ipynb` in the JupyterLab file browser.

## Dependencies

Managed in `pyproject.toml` and locked in `uv.lock`:

- **jupyterlab** — notebook interface
- **qiskit** — quantum circuit SDK
- **qiskit-aer** — local quantum simulator
- **cirq-core** — quantum circuit framework (Google)
- **stim** — fast stabilizer circuit simulator for QEC
- **pymatching** — decoder for the surface-code threshold experiment
- **matplotlib** — visualization
- **scipy** — scientific computing

To update dependencies:

```sh
uv lock --upgrade
```

## License

[The Unlicense](https://unlicense.org/) — public domain.
