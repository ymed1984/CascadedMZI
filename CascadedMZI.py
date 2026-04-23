from __future__ import annotations

import argparse
import csv
import math
from cmath import exp
from dataclasses import dataclass
from pathlib import Path


NM_PER_M = 1e9
UM_PER_M = 1e6
MIN_DB_FLOOR = -60.0


@dataclass(frozen=True)
class MZIStage:
    index: int
    fsr_nm: float
    delta_length_m: float
    coupler_ratio: float
    insertion_loss_db: float

    @property
    def delta_length_um(self) -> float:
        return self.delta_length_m * UM_PER_M


@dataclass(frozen=True)
class CascadedMZI:
    center_wavelength_nm: float
    effective_index: float
    group_index: float
    stages: tuple[MZIStage, ...]
    base_fsr_nm: float
    fsr_scale: float
    bend_radius_um: float
    extra_straight_um: float

    @property
    def total_delta_length_um(self) -> float:
        return sum(stage.delta_length_um for stage in self.stages)

    @property
    def estimated_footprint_um(self) -> float:
        # Very rough routing estimate:
        # each delay stage consumes two bends plus the excess straight section.
        bend_span = math.pi * self.bend_radius_um
        return sum(stage.delta_length_um + bend_span + self.extra_straight_um for stage in self.stages)


@dataclass(frozen=True)
class LatticeFilter:
    center_wavelength_nm: float
    effective_index: float
    group_index: float
    fsr_nm: float
    coupler_ratios: tuple[float, ...]
    delay_lengths_m: tuple[float, ...]
    preset_name: str

    @property
    def delay_lengths_um(self) -> tuple[float, ...]:
        return tuple(length * UM_PER_M for length in self.delay_lengths_m)

    @property
    def order(self) -> int:
        return len(self.delay_lengths_m)


def linear_to_db(power: float, floor_db: float = MIN_DB_FLOOR) -> float:
    if power <= 0.0:
        return floor_db
    return max(floor_db, 10.0 * math.log10(power))


def compute_delta_length_m(center_wavelength_nm: float, group_index: float, fsr_nm: float) -> float:
    wavelength_m = center_wavelength_nm / NM_PER_M
    fsr_m = fsr_nm / NM_PER_M
    return wavelength_m**2 / (group_index * fsr_m)


def compute_pi_length_m(center_wavelength_nm: float, effective_index: float) -> float:
    wavelength_m = center_wavelength_nm / NM_PER_M
    return wavelength_m / (2.0 * effective_index)


def get_lattice_filter_spec(
    preset_name: str,
    unit_delay_m: float,
    pi_length_m: float,
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    normalized_presets: dict[str, tuple[tuple[float, ...], tuple[float, ...]]] = {
        "dwivedi-flat": (
            (0.5, 0.13, 0.12, 0.5, 0.25),
            (unit_delay_m, 2.0 * unit_delay_m, -(2.0 * unit_delay_m + pi_length_m), -2.0 * unit_delay_m),
        ),
        "maxflat-n2": (
            tuple(math.sin(value * math.pi) ** 2 for value in (0.0833, 0.333, 0.25)),
            (2.0 * unit_delay_m - pi_length_m, unit_delay_m),
        ),
        "maxflat-n3": (
            tuple(math.sin(value * math.pi) ** 2 for value in (0.4664, 0.1591, 0.3755, 0.25)),
            (2.0 * unit_delay_m, 2.0 * unit_delay_m - pi_length_m, unit_delay_m),
        ),
        "maxflat-n4": (
            tuple(math.sin(value * math.pi) ** 2 for value in (0.3720, 0.2860, 0.4547, 0.1187, 0.25)),
            (2.0 * unit_delay_m - pi_length_m, 2.0 * unit_delay_m, 2.0 * unit_delay_m - pi_length_m, unit_delay_m),
        ),
        "cheby-n4": (
            tuple(math.sin(value * math.pi) ** 2 for value in (0.3498, 0.2448, 0.4186, 0.0797, 0.25)),
            (2.0 * unit_delay_m - pi_length_m, 2.0 * unit_delay_m, 2.0 * unit_delay_m - pi_length_m, unit_delay_m),
        ),
    }

    try:
        return normalized_presets[preset_name]
    except KeyError as exc:
        available = ", ".join(sorted(normalized_presets))
        raise ValueError(f"Unknown lattice preset '{preset_name}'. Available presets: {available}.") from exc


def build_lattice_filter(
    center_wavelength_nm: float,
    effective_index: float,
    group_index: float,
    fsr_nm: float,
    preset_name: str,
) -> LatticeFilter:
    if center_wavelength_nm <= 0 or fsr_nm <= 0:
        raise ValueError("Wavelength and FSR must be positive.")
    if group_index <= 0 or effective_index <= 0:
        raise ValueError("Indices must be positive.")

    unit_delay_m = compute_delta_length_m(center_wavelength_nm=center_wavelength_nm, group_index=group_index, fsr_nm=fsr_nm)
    pi_length_m = compute_pi_length_m(center_wavelength_nm=center_wavelength_nm, effective_index=effective_index)
    coupler_ratios, delay_lengths_m = get_lattice_filter_spec(
        preset_name=preset_name,
        unit_delay_m=unit_delay_m,
        pi_length_m=pi_length_m,
    )

    return LatticeFilter(
        center_wavelength_nm=center_wavelength_nm,
        effective_index=effective_index,
        group_index=group_index,
        fsr_nm=fsr_nm,
        coupler_ratios=coupler_ratios,
        delay_lengths_m=delay_lengths_m,
        preset_name=preset_name,
    )


def build_cascaded_mzi(
    center_wavelength_nm: float,
    effective_index: float,
    group_index: float,
    stage_count: int,
    base_fsr_nm: float,
    fsr_scale: float,
    coupler_ratio: float,
    insertion_loss_db: float,
    bend_radius_um: float,
    extra_straight_um: float,
) -> CascadedMZI:
    if stage_count < 1:
        raise ValueError("stage_count must be 1 or larger.")
    if center_wavelength_nm <= 0 or base_fsr_nm <= 0:
        raise ValueError("Wavelength and FSR must be positive.")
    if group_index <= 0 or effective_index <= 0:
        raise ValueError("Indices must be positive.")
    if fsr_scale <= 0:
        raise ValueError("fsr_scale must be positive.")
    if not 0 < coupler_ratio < 1:
        raise ValueError("coupler_ratio must be between 0 and 1.")

    stages: list[MZIStage] = []
    for stage_index in range(stage_count):
        fsr_nm = base_fsr_nm * (fsr_scale**stage_index)
        delta_length_m = compute_delta_length_m(center_wavelength_nm, group_index, fsr_nm)
        stages.append(
            MZIStage(
                index=stage_index + 1,
                fsr_nm=fsr_nm,
                delta_length_m=delta_length_m,
                coupler_ratio=coupler_ratio,
                insertion_loss_db=insertion_loss_db,
            )
        )

    return CascadedMZI(
        center_wavelength_nm=center_wavelength_nm,
        effective_index=effective_index,
        group_index=group_index,
        stages=tuple(stages),
        base_fsr_nm=base_fsr_nm,
        fsr_scale=fsr_scale,
        bend_radius_um=bend_radius_um,
        extra_straight_um=extra_straight_um,
    )


def stage_transfer(stage: MZIStage, wavelength_nm: float, center_wavelength_nm: float, effective_index: float) -> float:
    wavelength_m = wavelength_nm / NM_PER_M
    phase = 2.0 * math.pi * effective_index * stage.delta_length_m / wavelength_m

    kappa = stage.coupler_ratio
    # Simple balanced MZI intensity response with non-ideal coupler ratio.
    visibility = 4.0 * kappa * (1.0 - kappa)
    ideal_transmission = 1.0 - visibility * (math.sin(phase / 2.0) ** 2)
    loss_linear = 10.0 ** (-stage.insertion_loss_db / 10.0)

    # Anchor the phase around the requested center wavelength for readability.
    center_phase = 2.0 * math.pi * effective_index * stage.delta_length_m / (center_wavelength_nm / NM_PER_M)
    normalized = 0.5 + 0.5 * math.cos(phase - center_phase)
    return max(0.0, min(1.0, loss_linear * (0.5 * ideal_transmission + 0.5 * normalized)))


def coupler_matrix(coupling_ratio: float) -> tuple[tuple[complex, complex], tuple[complex, complex]]:
    through = math.sqrt(1.0 - coupling_ratio)
    cross = -1j * math.sqrt(coupling_ratio)
    return ((through, cross), (cross, through))


def multiply_matrix_vector(
    matrix: tuple[tuple[complex, complex], tuple[complex, complex]],
    vector: tuple[complex, complex],
) -> tuple[complex, complex]:
    return (
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1],
    )


def lattice_transfer(
    design: LatticeFilter,
    wavelength_nm: float,
) -> tuple[float, float]:
    wavelength_m = wavelength_nm / NM_PER_M
    field = (1.0 + 0.0j, 0.0 + 0.0j)

    for index, coupling_ratio in enumerate(design.coupler_ratios):
        field = multiply_matrix_vector(coupler_matrix(coupling_ratio), field)
        if index < len(design.delay_lengths_m):
            phase = 2.0 * math.pi * design.effective_index * design.delay_lengths_m[index] / wavelength_m
            field = (field[0] * exp(-1j * phase), field[1])

    through_power = abs(field[0]) ** 2
    cross_power = abs(field[1]) ** 2
    return through_power, cross_power


def spectrum(
    design: CascadedMZI,
    start_nm: float,
    stop_nm: float,
    points: int,
) -> list[tuple[float, float]]:
    if points < 2:
        raise ValueError("points must be 2 or larger.")
    if stop_nm <= start_nm:
        raise ValueError("stop_nm must be larger than start_nm.")

    values: list[tuple[float, float]] = []
    step = (stop_nm - start_nm) / (points - 1)
    for idx in range(points):
        wavelength_nm = start_nm + idx * step
        transmission = 1.0
        for stage in design.stages:
            transmission *= stage_transfer(
                stage=stage,
                wavelength_nm=wavelength_nm,
                center_wavelength_nm=design.center_wavelength_nm,
                effective_index=design.effective_index,
            )
        values.append((wavelength_nm, transmission))
    return values


def lattice_spectrum(
    design: LatticeFilter,
    start_nm: float,
    stop_nm: float,
    points: int,
) -> list[tuple[float, float, float]]:
    if points < 2:
        raise ValueError("points must be 2 or larger.")
    if stop_nm <= start_nm:
        raise ValueError("stop_nm must be larger than start_nm.")

    values: list[tuple[float, float, float]] = []
    step = (stop_nm - start_nm) / (points - 1)
    for idx in range(points):
        wavelength_nm = start_nm + idx * step
        through_power, cross_power = lattice_transfer(design, wavelength_nm)
        values.append((wavelength_nm, through_power, cross_power))
    return values


def save_spectrum_csv(path: Path, rows: list[tuple[float, float]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["wavelength_nm", "transmission_db"])
        writer.writerows((wavelength_nm, linear_to_db(transmission)) for wavelength_nm, transmission in rows)


def save_lattice_spectrum_csv(path: Path, rows: list[tuple[float, float, float]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["wavelength_nm", "through_db", "cross_db"])
        writer.writerows(
            (wavelength_nm, linear_to_db(through), linear_to_db(cross))
            for wavelength_nm, through, cross in rows
        )


def plot_spectrum(
    design: CascadedMZI,
    rows: list[tuple[float, float]],
    output_path: Path | None = None,
    show_plot: bool = False,
) -> None:
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "matplotlib is not installed. Install it with 'pip install matplotlib' to enable plotting."
        ) from exc

    wavelengths = [row[0] for row in rows]
    transmissions_db = [linear_to_db(row[1]) for row in rows]

    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot(wavelengths, transmissions_db, color="tab:blue", linewidth=2.0)
    ax.set_title(f"Cascaded MZI Spectrum ({len(design.stages)} stages)")
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("Transmission [dB]")
    ax.set_ylim(MIN_DB_FLOOR, 1.0)
    ax.grid(True, alpha=0.3)

    summary = (
        f"lambda0 = {design.center_wavelength_nm:.1f} nm\n"
        f"ng = {design.group_index:.2f}, neff = {design.effective_index:.2f}\n"
        f"base FSR = {design.base_fsr_nm:.2f} nm, scale = {design.fsr_scale:.2f}"
    )
    ax.text(
        0.02,
        0.98,
        summary,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.8"},
    )

    fig.tight_layout()

    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=200, bbox_inches="tight")

    if show_plot:
        plt.show()

    plt.close(fig)


def plot_lattice_spectrum(
    design: LatticeFilter,
    rows: list[tuple[float, float, float]],
    output_path: Path | None = None,
    show_plot: bool = False,
) -> None:
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "matplotlib is not installed. Install it with 'pip install matplotlib' to enable plotting."
        ) from exc

    wavelengths = [row[0] for row in rows]
    through_db = [linear_to_db(row[1]) for row in rows]
    cross_db = [linear_to_db(row[2]) for row in rows]

    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    ax.plot(wavelengths, through_db, label="Through", color="tab:blue", linewidth=2.0)
    ax.plot(wavelengths, cross_db, label="Cross", color="tab:orange", linewidth=2.0)
    ax.set_title(f"MZI Lattice Spectrum ({design.preset_name})")
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("Transmission [dB]")
    ax.set_ylim(MIN_DB_FLOOR, 1.0)
    ax.grid(True, alpha=0.3)
    ax.legend()

    summary = (
        f"lambda0 = {design.center_wavelength_nm:.1f} nm\n"
        f"FSR = {design.fsr_nm:.2f} nm\n"
        f"preset = {design.preset_name}\n"
        f"order = {design.order}"
    )
    ax.text(
        0.02,
        0.98,
        summary,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.8"},
    )

    fig.tight_layout()

    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=200, bbox_inches="tight")

    if show_plot:
        plt.show()

    plt.close(fig)


def format_report(design: CascadedMZI) -> str:
    lines = [
        "Cascaded MZI Simple Design Report",
        "=" * 33,
        f"Center wavelength : {design.center_wavelength_nm:.3f} nm",
        f"Effective index   : {design.effective_index:.4f}",
        f"Group index       : {design.group_index:.4f}",
        f"Stages            : {len(design.stages)}",
        f"Base FSR          : {design.base_fsr_nm:.3f} nm",
        f"FSR scale         : {design.fsr_scale:.3f}",
        f"Total dL          : {design.total_delta_length_um:.2f} um",
        f"Rough footprint   : {design.estimated_footprint_um:.2f} um",
        "",
        "Stage details",
        "-" * 33,
        f"{'Stage':>5} {'FSR [nm]':>10} {'dL [um]':>12} {'Coupler':>10} {'Loss [dB]':>10}",
    ]

    for stage in design.stages:
        lines.append(
            f"{stage.index:>5} "
            f"{stage.fsr_nm:>10.4f} "
            f"{stage.delta_length_um:>12.2f} "
            f"{stage.coupler_ratio:>10.2f} "
            f"{stage.insertion_loss_db:>10.2f}"
        )

    return "\n".join(lines)


def format_lattice_report(design: LatticeFilter) -> str:
    lines = [
        "MZI Lattice Filter Design Report",
        "=" * 32,
        f"Preset            : {design.preset_name}",
        f"Center wavelength : {design.center_wavelength_nm:.3f} nm",
        f"Effective index   : {design.effective_index:.4f}",
        f"Group index       : {design.group_index:.4f}",
        f"FSR               : {design.fsr_nm:.3f} nm",
        f"Order             : {design.order}",
        "",
        "Coupler details",
        "-" * 32,
        f"{'Index':>5} {'Power coupling':>16}",
    ]
    for index, coupling_ratio in enumerate(design.coupler_ratios, start=1):
        lines.append(f"{index:>5} {coupling_ratio:>16.6f}")

    lines.extend(
        [
            "",
            "Delay details",
            "-" * 32,
            f"{'Index':>5} {'dL [um]':>14}",
        ]
    )
    for index, delay_length_um in enumerate(design.delay_lengths_um, start=1):
        lines.append(f"{index:>5} {delay_length_um:>14.3f}")

    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Simple design helper for cascaded Mach-Zehnder interferometers."
    )
    parser.add_argument(
        "--architecture",
        choices=("simple", "lattice"),
        default="simple",
        help="Filter architecture to synthesize.",
    )
    parser.add_argument("--center-wavelength", type=float, default=1550.0, help="Center wavelength in nm.")
    parser.add_argument("--effective-index", type=float, default=2.40, help="Effective index.")
    parser.add_argument("--group-index", type=float, default=4.20, help="Group index.")
    parser.add_argument("--stages", type=int, default=3, help="Number of cascaded MZI stages.")
    parser.add_argument("--base-fsr", type=float, default=8.0, help="First-stage FSR in nm.")
    parser.add_argument(
        "--lattice-fsr",
        type=float,
        default=10.0,
        help="FSR in nm for lattice filter synthesis.",
    )
    parser.add_argument(
        "--lattice-preset",
        type=str,
        default="dwivedi-flat",
        help="Lattice preset: dwivedi-flat, maxflat-n2, maxflat-n3, maxflat-n4, cheby-n4.",
    )
    parser.add_argument(
        "--fsr-scale",
        type=float,
        default=2.0,
        help="Multiplier applied to each next-stage FSR.",
    )
    parser.add_argument(
        "--coupler-ratio",
        type=float,
        default=0.5,
        help="Power coupling ratio for each coupler (0 < ratio < 1).",
    )
    parser.add_argument(
        "--insertion-loss",
        type=float,
        default=0.2,
        help="Insertion loss per stage in dB.",
    )
    parser.add_argument("--bend-radius", type=float, default=30.0, help="Routing bend radius in um.")
    parser.add_argument(
        "--extra-straight",
        type=float,
        default=20.0,
        help="Extra straight routing margin per stage in um.",
    )
    parser.add_argument(
        "--span",
        type=float,
        default=20.0,
        help="Spectrum span around the center wavelength in nm.",
    )
    parser.add_argument(
        "--start-wavelength",
        type=float,
        default=None,
        help="Absolute spectrum start wavelength in nm. Overrides --span when set with --stop-wavelength.",
    )
    parser.add_argument(
        "--stop-wavelength",
        type=float,
        default=None,
        help="Absolute spectrum stop wavelength in nm. Overrides --span when set with --start-wavelength.",
    )
    parser.add_argument("--points", type=int, default=401, help="Number of sampled spectrum points.")
    parser.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="Optional output path for sampled spectrum CSV.",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Display the spectrum using matplotlib.",
    )
    parser.add_argument(
        "--plot-file",
        type=Path,
        default=None,
        help="Optional output path for a plotted spectrum image such as spectrum.png.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    has_start = args.start_wavelength is not None
    has_stop = args.stop_wavelength is not None
    if has_start != has_stop:
        parser.exit(status=2, message="Both --start-wavelength and --stop-wavelength must be provided together.\n")

    if has_start and has_stop:
        start_nm = args.start_wavelength
        stop_nm = args.stop_wavelength
    else:
        start_nm = args.center_wavelength - args.span / 2.0
        stop_nm = args.center_wavelength + args.span / 2.0

    if args.architecture == "simple":
        design = build_cascaded_mzi(
            center_wavelength_nm=args.center_wavelength,
            effective_index=args.effective_index,
            group_index=args.group_index,
            stage_count=args.stages,
            base_fsr_nm=args.base_fsr,
            fsr_scale=args.fsr_scale,
            coupler_ratio=args.coupler_ratio,
            insertion_loss_db=args.insertion_loss,
            bend_radius_um=args.bend_radius,
            extra_straight_um=args.extra_straight,
        )

        print(format_report(design))
        rows = spectrum(design, start_nm=start_nm, stop_nm=stop_nm, points=args.points)
        peak = max(rows, key=lambda item: item[1])
        valley = min(rows, key=lambda item: item[1])
        print()
        print("Spectrum summary")
        print("-" * 33)
        print(f"Sample range      : {start_nm:.3f} - {stop_nm:.3f} nm")
        print(f"Peak transmission : {linear_to_db(peak[1]):.2f} dB @ {peak[0]:.3f} nm")
        print(f"Valley            : {linear_to_db(valley[1]):.2f} dB @ {valley[0]:.3f} nm")

        if args.csv is not None:
            args.csv.parent.mkdir(parents=True, exist_ok=True)
            save_spectrum_csv(args.csv, rows)
            print(f"CSV saved         : {args.csv}")

        if args.plot or args.plot_file is not None:
            try:
                plot_spectrum(
                    design=design,
                    rows=rows,
                    output_path=args.plot_file,
                    show_plot=args.plot,
                )
            except RuntimeError as exc:
                parser.exit(status=1, message=f"{exc}\n")
            if args.plot_file is not None:
                print(f"Plot saved        : {args.plot_file}")
        return

    design = build_lattice_filter(
        center_wavelength_nm=args.center_wavelength,
        effective_index=args.effective_index,
        group_index=args.group_index,
        fsr_nm=args.lattice_fsr,
        preset_name=args.lattice_preset,
    )
    print(format_lattice_report(design))
    rows = lattice_spectrum(design, start_nm=start_nm, stop_nm=stop_nm, points=args.points)
    peak_through = max(rows, key=lambda item: item[1])
    peak_cross = max(rows, key=lambda item: item[2])
    min_cross = min(rows, key=lambda item: item[2])
    print()
    print("Spectrum summary")
    print("-" * 32)
    print(f"Sample range      : {start_nm:.3f} - {stop_nm:.3f} nm")
    print(f"Peak through      : {linear_to_db(peak_through[1]):.2f} dB @ {peak_through[0]:.3f} nm")
    print(f"Peak cross        : {linear_to_db(peak_cross[2]):.2f} dB @ {peak_cross[0]:.3f} nm")
    print(f"Min cross         : {linear_to_db(min_cross[2]):.2f} dB @ {min_cross[0]:.3f} nm")

    if args.csv is not None:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        save_lattice_spectrum_csv(args.csv, rows)
        print(f"CSV saved         : {args.csv}")

    if args.plot or args.plot_file is not None:
        try:
            plot_lattice_spectrum(
                design=design,
                rows=rows,
                output_path=args.plot_file,
                show_plot=args.plot,
            )
        except RuntimeError as exc:
            parser.exit(status=1, message=f"{exc}\n")
        if args.plot_file is not None:
            print(f"Plot saved        : {args.plot_file}")


if __name__ == "__main__":
    main()
