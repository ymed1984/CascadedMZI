from __future__ import annotations

import argparse
import bisect
import csv
import math
from cmath import exp
from dataclasses import dataclass
from pathlib import Path


NM_PER_M = 1e9
UM_PER_M = 1e6
C_M_PER_S = 299_792_458.0
MIN_DB_FLOOR = -60.0
DEFAULT_CENTER_WAVELENGTH_NM = 1550.0
DEFAULT_EFFECTIVE_INDEX = 2.40
DEFAULT_GROUP_INDEX = 4.20
PAPER_LAN_CENTER_WAVELENGTH_NM = 1291.0
PAPER_LAN_EFFECTIVE_INDEX = 2.93
PAPER_LAN_GROUP_INDEX = 3.83
PAPER_LAN_CHANNEL_WAVELENGTHS_NM = (1273.5, 1277.9, 1282.3, 1286.7, 1295.6, 1300.1, 1304.6, 1309.1)
PUB_16WDM_CHANNEL_COUNT = 16
PUB_16WDM_CHANNEL_SPACING_GHZ = 800.0
PUB_16WDM_COUPLER_EXCESS_LOSS_DB = 0.03
PUB_16WDM_PORT_CHANNELS = (1, 9, 13, 5, 3, 11, 15, 7, 2, 10, 14, 6, 4, 12, 16, 8)


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


@dataclass(frozen=True)
class DemuxMZIStage:
    name: str
    delta_length_m: float
    coupler_ratio: float
    phase_offset_rad: float

    @property
    def delta_length_um(self) -> float:
        return self.delta_length_m * UM_PER_M


@dataclass(frozen=True)
class LANDemux:
    center_wavelength_nm: float
    effective_index: float
    group_index: float
    channel_spacing_nm: float
    stages: tuple[DemuxMZIStage, ...]
    channel_wavelengths_nm: tuple[float, ...]
    port_wavelengths_nm: tuple[float, ...]

    @property
    def order(self) -> int:
        return 3


@dataclass(frozen=True)
class CouplerResponse:
    wavelengths_nm: tuple[float, ...]
    coupling_ratios: tuple[float, ...]
    excess_loss_db: tuple[float, ...]
    source_name: str = "constant"

    def coupling_ratio(self, wavelength_nm: float) -> float:
        return interpolate_series(self.wavelengths_nm, self.coupling_ratios, wavelength_nm)

    def loss_db(self, wavelength_nm: float) -> float:
        return interpolate_series(self.wavelengths_nm, self.excess_loss_db, wavelength_nm)


@dataclass(frozen=True)
class WDMSplitterStage:
    name: str
    level: int
    node_index: int
    delta_length_m: float
    phase_offset_rad: float
    through_channels: tuple[int, ...]
    cross_channels: tuple[int, ...]

    @property
    def delta_length_um(self) -> float:
        return self.delta_length_m * UM_PER_M


@dataclass(frozen=True)
class CascadedWDMDemux:
    center_wavelength_nm: float
    effective_index: float
    group_index: float
    channel_spacing_ghz: float
    channel_wavelengths_nm: tuple[float, ...]
    port_channels: tuple[int, ...]
    stages: tuple[WDMSplitterStage, ...]
    coupler_response: CouplerResponse
    source_name: str

    @property
    def channel_count(self) -> int:
        return len(self.channel_wavelengths_nm)

    @property
    def order(self) -> int:
        return int(math.log2(self.channel_count))


def linear_to_db(power: float, floor_db: float = MIN_DB_FLOOR) -> float:
    if power <= 0.0:
        return floor_db
    return max(floor_db, 10.0 * math.log10(power))


def interpolate_series(x_values: tuple[float, ...], y_values: tuple[float, ...], x: float) -> float:
    if len(x_values) != len(y_values) or not x_values:
        raise ValueError("Interpolation series must have the same nonzero length.")
    if len(x_values) == 1 or x <= x_values[0]:
        return y_values[0]
    if x >= x_values[-1]:
        return y_values[-1]

    upper_index = bisect.bisect_left(x_values, x)
    x0 = x_values[upper_index - 1]
    x1 = x_values[upper_index]
    y0 = y_values[upper_index - 1]
    y1 = y_values[upper_index]
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)


def constant_coupler_response(
    coupling_ratio: float,
    excess_loss_db: float = 0.0,
    source_name: str = "constant",
) -> CouplerResponse:
    if not 0 < coupling_ratio < 1:
        raise ValueError("coupling_ratio must be between 0 and 1.")
    if excess_loss_db < 0:
        raise ValueError("excess_loss_db must not be negative.")
    return CouplerResponse(
        wavelengths_nm=(0.0,),
        coupling_ratios=(coupling_ratio,),
        excess_loss_db=(excess_loss_db,),
        source_name=source_name,
    )


def load_coupler_response_csv(
    path: Path,
    ratio_column: str,
    wavelength_column: str,
    loss_column: str | None,
    fallback_loss_db: float,
) -> CouplerResponse:
    rows: list[tuple[float, float, float]] = []
    with path.open(newline="", encoding="utf-8-sig") as csv_file:
        reader = csv.DictReader(csv_file)
        if reader.fieldnames is None:
            raise ValueError(f"Coupler file '{path}' has no header row.")
        missing = [column for column in (wavelength_column, ratio_column) if column not in reader.fieldnames]
        if missing:
            raise ValueError(f"Coupler file '{path}' is missing columns: {', '.join(missing)}.")
        has_loss_column = loss_column is not None and loss_column in reader.fieldnames
        for row in reader:
            wavelength_nm = float(row[wavelength_column])
            coupling_ratio = float(row[ratio_column])
            loss_db = float(row[loss_column]) if has_loss_column else fallback_loss_db
            if not 0 < coupling_ratio < 1:
                raise ValueError(f"Coupler ratio must be between 0 and 1, got {coupling_ratio}.")
            if loss_db < 0:
                raise ValueError(f"Coupler excess loss must not be negative, got {loss_db}.")
            rows.append((wavelength_nm, coupling_ratio, loss_db))

    if not rows:
        raise ValueError(f"Coupler file '{path}' does not contain data rows.")

    rows.sort(key=lambda item: item[0])
    return CouplerResponse(
        wavelengths_nm=tuple(row[0] for row in rows),
        coupling_ratios=tuple(row[1] for row in rows),
        excess_loss_db=tuple(row[2] for row in rows),
        source_name=str(path),
    )


def compute_delta_length_m(center_wavelength_nm: float, group_index: float, fsr_nm: float) -> float:
    wavelength_m = center_wavelength_nm / NM_PER_M
    fsr_m = fsr_nm / NM_PER_M
    return wavelength_m**2 / (group_index * fsr_m)


def compute_delta_length_from_fsr_hz(group_index: float, fsr_hz: float) -> float:
    return C_M_PER_S / (group_index * fsr_hz)


def compute_pi_length_m(center_wavelength_nm: float, effective_index: float) -> float:
    wavelength_m = center_wavelength_nm / NM_PER_M
    return wavelength_m / (2.0 * effective_index)


def amzi_power_from_phase(coupling_ratio: float, phase_rad: float) -> tuple[float, float]:
    field = (1.0 + 0.0j, 0.0 + 0.0j)
    coupler = coupler_matrix(coupling_ratio)
    field = multiply_matrix_vector(coupler, field)
    field = (field[0] * exp(-1j * phase_rad), field[1])
    field = multiply_matrix_vector(coupler, field)
    return abs(field[0]) ** 2, abs(field[1]) ** 2


def amzi_phase(
    delta_length_m: float,
    wavelength_nm: float,
    center_wavelength_nm: float,
    effective_index: float,
    group_index: float,
    phase_offset_rad: float | None = None,
) -> float:
    wavelength_m = wavelength_nm / NM_PER_M
    center_wavelength_m = center_wavelength_nm / NM_PER_M
    if phase_offset_rad is None:
        phase = 2.0 * math.pi * effective_index * delta_length_m / center_wavelength_m
    else:
        phase = phase_offset_rad
    return phase + 2.0 * math.pi * group_index * delta_length_m * (1.0 / wavelength_m - 1.0 / center_wavelength_m)


def optimize_amzi_phase_offset(
    delta_length_m: float,
    coupling_ratio: float,
    through_wavelengths_nm: tuple[float, ...],
    cross_wavelengths_nm: tuple[float, ...],
    center_wavelength_nm: float,
    effective_index: float,
    group_index: float,
) -> float:
    best_offset = 0.0
    best_score = -math.inf
    for step in range(720):
        offset = 2.0 * math.pi * step / 720
        score = 0.0
        for wavelength_nm in through_wavelengths_nm:
            phase = amzi_phase(
                delta_length_m,
                wavelength_nm,
                center_wavelength_nm,
                effective_index,
                group_index,
                offset,
            )
            through_power, _ = amzi_power_from_phase(coupling_ratio, phase)
            score += math.log(max(through_power, 1e-12))
        for wavelength_nm in cross_wavelengths_nm:
            phase = amzi_phase(
                delta_length_m,
                wavelength_nm,
                center_wavelength_nm,
                effective_index,
                group_index,
                offset,
            )
            _, cross_power = amzi_power_from_phase(coupling_ratio, phase)
            score += math.log(max(cross_power, 1e-12))
        if score > best_score:
            best_score = score
            best_offset = offset
    return best_offset


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


def build_lan_wdm_demux(
    center_wavelength_nm: float,
    effective_index: float,
    group_index: float,
    channel_spacing_nm: float,
    coupler_ratio: float,
    channel_wavelengths_nm: tuple[float, ...] = PAPER_LAN_CHANNEL_WAVELENGTHS_NM,
    optimize_phase_offsets: bool = True,
) -> LANDemux:
    if center_wavelength_nm <= 0 or channel_spacing_nm <= 0:
        raise ValueError("Wavelength and channel spacing must be positive.")
    if group_index <= 0 or effective_index <= 0:
        raise ValueError("Indices must be positive.")
    if not 0 < coupler_ratio < 1:
        raise ValueError("coupler_ratio must be between 0 and 1.")
    if len(channel_wavelengths_nm) != 8:
        raise ValueError("LAN-WDM demux requires exactly eight target channel wavelengths.")

    # Photonics 2022, 9, 252 uses FSR1 = 2 * channel spacing for the first splitter.
    delta_l1_m = compute_delta_length_m(
        center_wavelength_nm=center_wavelength_nm,
        group_index=group_index,
        fsr_nm=2.0 * channel_spacing_nm,
    )
    delta_l_shift_m = compute_pi_length_m(center_wavelength_nm, effective_index) * 2.0

    channel_1, channel_2, channel_3, channel_4, channel_6, channel_7, channel_8, channel_9 = channel_wavelengths_nm
    stage_specs = (
        ("1st", delta_l1_m, (channel_1, channel_3, channel_7, channel_9), (channel_2, channel_4, channel_6, channel_8)),
        ("2nd_A", delta_l1_m / 2.0, (channel_1, channel_9), (channel_3, channel_7)),
        (
            "2nd_B",
            delta_l1_m / 2.0 + 0.75 * delta_l_shift_m,
            (channel_2, channel_6),
            (channel_4, channel_8),
        ),
        ("3rd_A", delta_l1_m / 8.0, (channel_1,), (channel_9,)),
        (
            "3rd_B",
            delta_l1_m / 4.0 + 0.25 * delta_l_shift_m,
            (channel_3,),
            (channel_7,),
        ),
        (
            "3rd_C",
            delta_l1_m / 4.0 + 0.125 * delta_l_shift_m,
            (channel_2,),
            (channel_6,),
        ),
        (
            "3rd_D",
            delta_l1_m / 4.0 + 0.375 * delta_l_shift_m,
            (channel_4,),
            (channel_8,),
        ),
    )

    return LANDemux(
        center_wavelength_nm=center_wavelength_nm,
        effective_index=effective_index,
        group_index=group_index,
        channel_spacing_nm=channel_spacing_nm,
        stages=tuple(
            DemuxMZIStage(
                name=name,
                delta_length_m=delta_length_m,
                coupler_ratio=coupler_ratio,
                phase_offset_rad=(
                    optimize_amzi_phase_offset(
                        delta_length_m=delta_length_m,
                        coupling_ratio=coupler_ratio,
                        through_wavelengths_nm=tuple(through_wavelengths_nm),
                        cross_wavelengths_nm=tuple(cross_wavelengths_nm),
                        center_wavelength_nm=center_wavelength_nm,
                        effective_index=effective_index,
                        group_index=group_index,
                    )
                    if optimize_phase_offsets
                    else 0.0
                ),
            )
            for name, delta_length_m, through_wavelengths_nm, cross_wavelengths_nm in stage_specs
        ),
        channel_wavelengths_nm=channel_wavelengths_nm,
        port_wavelengths_nm=(channel_1, channel_9, channel_3, channel_7, channel_2, channel_6, channel_4, channel_8),
    )


def frequency_grid_wavelengths_nm(
    center_wavelength_nm: float,
    channel_spacing_ghz: float,
    channel_count: int,
) -> tuple[float, ...]:
    if channel_count < 2 or channel_count & (channel_count - 1):
        raise ValueError("channel_count must be a power of two and at least 2.")
    if center_wavelength_nm <= 0 or channel_spacing_ghz <= 0:
        raise ValueError("Center wavelength and channel spacing must be positive.")

    center_frequency_hz = C_M_PER_S / (center_wavelength_nm / NM_PER_M)
    spacing_hz = channel_spacing_ghz * 1e9
    midpoint = (channel_count + 1) / 2.0
    wavelengths: list[float] = []
    for channel_index in range(1, channel_count + 1):
        frequency_hz = center_frequency_hz + (midpoint - channel_index) * spacing_hz
        wavelengths.append(C_M_PER_S / frequency_hz * NM_PER_M)
    return tuple(wavelengths)


def build_pub_16wdm_demux(
    center_wavelength_nm: float,
    effective_index: float,
    group_index: float,
    channel_spacing_ghz: float,
    coupler_response: CouplerResponse,
    optimize_phase_offsets: bool = True,
) -> CascadedWDMDemux:
    if center_wavelength_nm <= 0 or channel_spacing_ghz <= 0:
        raise ValueError("Center wavelength and channel spacing must be positive.")
    if effective_index <= 0 or group_index <= 0:
        raise ValueError("Indices must be positive.")

    channel_wavelengths_nm = frequency_grid_wavelengths_nm(
        center_wavelength_nm=center_wavelength_nm,
        channel_spacing_ghz=channel_spacing_ghz,
        channel_count=PUB_16WDM_CHANNEL_COUNT,
    )
    channel_wavelength_by_index = {
        channel_index: channel_wavelengths_nm[channel_index - 1]
        for channel_index in range(1, PUB_16WDM_CHANNEL_COUNT + 1)
    }
    base_delta_length_m = compute_delta_length_from_fsr_hz(
        group_index=group_index,
        fsr_hz=2.0 * channel_spacing_ghz * 1e9,
    )

    stages: list[WDMSplitterStage] = []

    def add_stages(channels: tuple[int, ...], level: int, node_index: int) -> None:
        if len(channels) == 1:
            return
        half = len(channels) // 2
        through_channels = channels[:half]
        cross_channels = channels[half:]
        delta_length_m = base_delta_length_m / (2 ** (level - 1))
        coupling_ratio = coupler_response.coupling_ratio(center_wavelength_nm)
        phase_offset_rad = (
            optimize_amzi_phase_offset(
                delta_length_m=delta_length_m,
                coupling_ratio=coupling_ratio,
                through_wavelengths_nm=tuple(channel_wavelength_by_index[index] for index in through_channels),
                cross_wavelengths_nm=tuple(channel_wavelength_by_index[index] for index in cross_channels),
                center_wavelength_nm=center_wavelength_nm,
                effective_index=effective_index,
                group_index=group_index,
            )
            if optimize_phase_offsets
            else 0.0
        )
        name = "MZI 1" if level == 1 else f"MZI {level}.{node_index}"
        stages.append(
            WDMSplitterStage(
                name=name,
                level=level,
                node_index=node_index,
                delta_length_m=delta_length_m,
                phase_offset_rad=phase_offset_rad,
                through_channels=through_channels,
                cross_channels=cross_channels,
            )
        )
        add_stages(through_channels, level + 1, 2 * node_index - 1)
        add_stages(cross_channels, level + 1, 2 * node_index)

    add_stages(PUB_16WDM_PORT_CHANNELS, level=1, node_index=1)
    stages.sort(key=lambda stage: (stage.level, stage.node_index))

    return CascadedWDMDemux(
        center_wavelength_nm=center_wavelength_nm,
        effective_index=effective_index,
        group_index=group_index,
        channel_spacing_ghz=channel_spacing_ghz,
        channel_wavelengths_nm=channel_wavelengths_nm,
        port_channels=PUB_16WDM_PORT_CHANNELS,
        stages=tuple(stages),
        coupler_response=coupler_response,
        source_name="pub_5335 Fig. 5 provisional 16-channel WDM",
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


def amzi_port_powers(
    delta_length_m: float,
    coupling_ratio: float,
    wavelength_nm: float,
    center_wavelength_nm: float,
    effective_index: float,
    group_index: float,
    center_phase_rad: float | None = None,
) -> tuple[float, float]:
    phase = amzi_phase(
        delta_length_m=delta_length_m,
        wavelength_nm=wavelength_nm,
        center_wavelength_nm=center_wavelength_nm,
        effective_index=effective_index,
        group_index=group_index,
        phase_offset_rad=center_phase_rad,
    )
    return amzi_power_from_phase(coupling_ratio, phase)


def amzi_port_powers_with_coupler(
    delta_length_m: float,
    coupler_response: CouplerResponse,
    wavelength_nm: float,
    center_wavelength_nm: float,
    effective_index: float,
    group_index: float,
    center_phase_rad: float,
) -> tuple[float, float]:
    through_power, cross_power = amzi_port_powers(
        delta_length_m=delta_length_m,
        coupling_ratio=coupler_response.coupling_ratio(wavelength_nm),
        wavelength_nm=wavelength_nm,
        center_wavelength_nm=center_wavelength_nm,
        effective_index=effective_index,
        group_index=group_index,
        center_phase_rad=center_phase_rad,
    )
    mzi_loss_linear = 10.0 ** (-(2.0 * coupler_response.loss_db(wavelength_nm)) / 10.0)
    return mzi_loss_linear * through_power, mzi_loss_linear * cross_power


def stage_transfer(
    stage: MZIStage,
    wavelength_nm: float,
    center_wavelength_nm: float,
    effective_index: float,
    group_index: float,
) -> float:
    # Use the constructive AMZI output port as the stage transmission.
    loss_linear = 10.0 ** (-stage.insertion_loss_db / 10.0)
    _, cross_power = amzi_port_powers(
        delta_length_m=stage.delta_length_m,
        coupling_ratio=stage.coupler_ratio,
        wavelength_nm=wavelength_nm,
        center_wavelength_nm=center_wavelength_nm,
        effective_index=effective_index,
        group_index=group_index,
        center_phase_rad=0.0,
    )
    return max(0.0, min(1.0, loss_linear * cross_power))


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
                group_index=design.group_index,
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


def lan_demux_transfer(design: LANDemux, wavelength_nm: float) -> tuple[float, ...]:
    stages = {stage.name: stage for stage in design.stages}

    def split(stage_name: str, input_power: float) -> tuple[float, float]:
        stage = stages[stage_name]
        through_power, cross_power = amzi_port_powers(
            delta_length_m=stage.delta_length_m,
            coupling_ratio=stage.coupler_ratio,
            wavelength_nm=wavelength_nm,
            center_wavelength_nm=design.center_wavelength_nm,
            effective_index=design.effective_index,
            group_index=design.group_index,
            center_phase_rad=stage.phase_offset_rad,
        )
        return input_power * through_power, input_power * cross_power

    root_through, root_cross = split("1st", 1.0)
    stage_2a_through, stage_2a_cross = split("2nd_A", root_through)
    stage_2b_through, stage_2b_cross = split("2nd_B", root_cross)

    port_1, port_2 = split("3rd_A", stage_2a_through)
    port_3, port_4 = split("3rd_B", stage_2a_cross)
    port_5, port_6 = split("3rd_C", stage_2b_through)
    port_7, port_8 = split("3rd_D", stage_2b_cross)
    return (port_1, port_2, port_3, port_4, port_5, port_6, port_7, port_8)


def lan_demux_spectrum(
    design: LANDemux,
    start_nm: float,
    stop_nm: float,
    points: int,
) -> list[tuple[float, tuple[float, ...]]]:
    if points < 2:
        raise ValueError("points must be 2 or larger.")
    if stop_nm <= start_nm:
        raise ValueError("stop_nm must be larger than start_nm.")

    values: list[tuple[float, tuple[float, ...]]] = []
    step = (stop_nm - start_nm) / (points - 1)
    for idx in range(points):
        wavelength_nm = start_nm + idx * step
        values.append((wavelength_nm, lan_demux_transfer(design, wavelength_nm)))
    return values


def cascaded_wdm_transfer(design: CascadedWDMDemux, wavelength_nm: float) -> tuple[float, ...]:
    node_inputs: dict[tuple[int, int], float] = {(1, 1): 1.0}
    channel_powers: dict[int, float] = {}

    for stage in design.stages:
        input_power = node_inputs.get((stage.level, stage.node_index), 0.0)
        through_power, cross_power = amzi_port_powers_with_coupler(
            delta_length_m=stage.delta_length_m,
            coupler_response=design.coupler_response,
            wavelength_nm=wavelength_nm,
            center_wavelength_nm=design.center_wavelength_nm,
            effective_index=design.effective_index,
            group_index=design.group_index,
            center_phase_rad=stage.phase_offset_rad,
        )
        through_output = input_power * through_power
        cross_output = input_power * cross_power
        if stage.level < design.order:
            node_inputs[(stage.level + 1, 2 * stage.node_index - 1)] = through_output
            node_inputs[(stage.level + 1, 2 * stage.node_index)] = cross_output
        else:
            channel_powers[stage.through_channels[0]] = through_output
            channel_powers[stage.cross_channels[0]] = cross_output

    return tuple(channel_powers[channel] for channel in design.port_channels)


def cascaded_wdm_spectrum(
    design: CascadedWDMDemux,
    start_nm: float,
    stop_nm: float,
    points: int,
) -> list[tuple[float, tuple[float, ...]]]:
    if points < 2:
        raise ValueError("points must be 2 or larger.")
    if stop_nm <= start_nm:
        raise ValueError("stop_nm must be larger than start_nm.")

    values: list[tuple[float, tuple[float, ...]]] = []
    step = (stop_nm - start_nm) / (points - 1)
    for idx in range(points):
        wavelength_nm = start_nm + idx * step
        values.append((wavelength_nm, cascaded_wdm_transfer(design, wavelength_nm)))
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


def save_lan_demux_spectrum_csv(path: Path, rows: list[tuple[float, tuple[float, ...]]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["wavelength_nm", *(f"port_{index}_db" for index in range(1, 9))])
        writer.writerows(
            (wavelength_nm, *(linear_to_db(power) for power in port_powers))
            for wavelength_nm, port_powers in rows
        )


def save_multiport_spectrum_csv(path: Path, rows: list[tuple[float, tuple[float, ...]]], port_count: int) -> None:
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["wavelength_nm", *(f"port_{index}_db" for index in range(1, port_count + 1))])
        writer.writerows(
            (wavelength_nm, *(linear_to_db(power) for power in port_powers))
            for wavelength_nm, port_powers in rows
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
    ax.set_title(f"Cascaded AMZI Spectrum ({len(design.stages)} stages)")
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


def plot_lan_demux_spectrum(
    design: LANDemux,
    rows: list[tuple[float, tuple[float, ...]]],
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

    fig, ax = plt.subplots(figsize=(9.0, 5.2))
    for port_index in range(8):
        port_db = [linear_to_db(row[1][port_index]) for row in rows]
        ax.plot(wavelengths, port_db, label=f"Port {port_index + 1}", linewidth=1.7)

    ax.set_title("Three-Stage Cascaded MZI LAN-WDM Demux")
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("Transmission [dB]")
    ax.set_ylim(MIN_DB_FLOOR, 1.0)
    ax.grid(True, alpha=0.3)
    ax.legend(ncol=4, fontsize=8)

    summary = (
        f"lambda0 = {design.center_wavelength_nm:.1f} nm\n"
        f"spacing = {design.channel_spacing_nm:.2f} nm\n"
        f"ng = {design.group_index:.2f}, neff = {design.effective_index:.2f}"
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


def plot_cascaded_wdm_spectrum(
    design: CascadedWDMDemux,
    rows: list[tuple[float, tuple[float, ...]]],
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

    fig, ax = plt.subplots(figsize=(10.0, 5.6))
    for port_index, channel_index in enumerate(design.port_channels):
        port_db = [linear_to_db(row[1][port_index]) for row in rows]
        ax.plot(wavelengths, port_db, label=f"P{port_index + 1}/ch{channel_index}", linewidth=1.25)

    ax.set_title("16-Channel Cascaded MZI WDM Demux")
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("Transmission [dB]")
    ax.set_ylim(MIN_DB_FLOOR, 1.0)
    ax.grid(True, alpha=0.3)
    ax.legend(ncol=4, fontsize=7)

    summary = (
        f"lambda0 = {design.center_wavelength_nm:.1f} nm\n"
        f"spacing = {design.channel_spacing_ghz:.1f} GHz\n"
        f"coupler = {design.coupler_response.source_name}"
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
    title = "Cascaded AMZI Simple Design Report"
    lines = [
        title,
        "=" * len(title),
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


def format_lan_demux_report(design: LANDemux) -> str:
    title = "Three-Stage Cascaded MZI LAN-WDM Demux Report"
    lines = [
        title,
        "=" * len(title),
        f"Center wavelength : {design.center_wavelength_nm:.3f} nm",
        f"Effective index   : {design.effective_index:.4f}",
        f"Group index       : {design.group_index:.4f}",
        f"Channel spacing   : {design.channel_spacing_nm:.3f} nm",
        f"Target channels   : {', '.join(f'{wavelength_nm:.1f}' for wavelength_nm in design.channel_wavelengths_nm)} nm",
        f"Port targets      : {', '.join(f'P{index + 1}={wavelength_nm:.1f}' for index, wavelength_nm in enumerate(design.port_wavelengths_nm))} nm",
        f"Tree order        : {design.order}",
        "",
        "MZI delay settings",
        "-" * 45,
        f"{'MZI':>7} {'dL [um]':>14} {'Coupler':>10} {'Phase [rad]':>12}",
    ]
    for stage in design.stages:
        lines.append(
            f"{stage.name:>7} {stage.delta_length_um:>14.3f} "
            f"{stage.coupler_ratio:>10.2f} {stage.phase_offset_rad:>12.3f}"
        )
    return "\n".join(lines)


def format_cascaded_wdm_report(design: CascadedWDMDemux) -> str:
    title = "16-Channel Cascaded MZI WDM Demux Report"
    lines = [
        title,
        "=" * len(title),
        f"Source            : {design.source_name}",
        f"Center wavelength : {design.center_wavelength_nm:.3f} nm",
        f"Effective index   : {design.effective_index:.4f}",
        f"Group index       : {design.group_index:.4f}",
        f"Channel spacing   : {design.channel_spacing_ghz:.3f} GHz",
        f"Channels          : {design.channel_count}",
        f"Tree order        : {design.order}",
        f"Coupler source    : {design.coupler_response.source_name}",
        "",
        "Channel wavelengths",
        "-" * 64,
        f"{'Channel':>8} {'Wavelength [nm]':>18} {'Output port':>12}",
    ]
    port_by_channel = {channel: port_index for port_index, channel in enumerate(design.port_channels, start=1)}
    for channel_index, wavelength_nm in enumerate(design.channel_wavelengths_nm, start=1):
        lines.append(f"{channel_index:>8} {wavelength_nm:>18.3f} {port_by_channel[channel_index]:>12}")

    lines.extend(
        [
            "",
            "MZI delay settings",
            "-" * 80,
            f"{'MZI':>8} {'Level':>5} {'Node':>5} {'dL [um]':>12} {'Phase [rad]':>12} {'Through ch.':>18} {'Cross ch.':>18}",
        ]
    )
    for stage in design.stages:
        through = ",".join(str(channel) for channel in stage.through_channels)
        cross = ",".join(str(channel) for channel in stage.cross_channels)
        lines.append(
            f"{stage.name:>8} {stage.level:>5} {stage.node_index:>5} "
            f"{stage.delta_length_um:>12.3f} {stage.phase_offset_rad:>12.3f} "
            f"{through:>18} {cross:>18}"
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
        description="Simple design helper for cascaded asymmetric Mach-Zehnder interferometers."
    )
    parser.add_argument(
        "--architecture",
        choices=("simple", "lattice", "lan-wdm", "pub-16wdm"),
        default="simple",
        help="Filter architecture to synthesize.",
    )
    parser.add_argument("--center-wavelength", type=float, default=None, help="Center wavelength in nm.")
    parser.add_argument("--effective-index", type=float, default=None, help="Effective index.")
    parser.add_argument("--group-index", type=float, default=None, help="Group index.")
    parser.add_argument("--stages", type=int, default=3, help="Number of cascaded MZI stages.")
    parser.add_argument("--base-fsr", type=float, default=8.0, help="First-stage FSR in nm.")
    parser.add_argument(
        "--channel-spacing",
        type=float,
        default=4.5,
        help="LAN-WDM channel spacing in nm for the three-stage binary-tree demux.",
    )
    parser.add_argument(
        "--channel-spacing-ghz",
        type=float,
        default=PUB_16WDM_CHANNEL_SPACING_GHZ,
        help="Frequency-domain channel spacing in GHz for pub-16wdm.",
    )
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
        "--coupler-excess-loss",
        type=float,
        default=None,
        help="Per-coupler excess loss in dB for pub-16wdm. Defaults to 0.03 dB for pub-16wdm.",
    )
    parser.add_argument(
        "--coupler-file",
        type=Path,
        default=None,
        help="Optional Lumerical-style coupler CSV with wavelength_nm, coupling_ratio, and optional excess_loss_db.",
    )
    parser.add_argument(
        "--coupler-wavelength-column",
        type=str,
        default="wavelength_nm",
        help="Wavelength column name in --coupler-file.",
    )
    parser.add_argument(
        "--coupler-ratio-column",
        type=str,
        default="coupling_ratio",
        help="Power coupling ratio column name in --coupler-file.",
    )
    parser.add_argument(
        "--coupler-loss-column",
        type=str,
        default="excess_loss_db",
        help="Optional per-coupler excess loss column name in --coupler-file.",
    )
    parser.add_argument(
        "--no-phase-offset",
        action="store_true",
        help="Disable automatic phase-offset alignment for binary-tree demux architectures.",
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
        default=None,
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

    if args.architecture == "lan-wdm":
        center_wavelength_nm = args.center_wavelength or PAPER_LAN_CENTER_WAVELENGTH_NM
        effective_index = args.effective_index or PAPER_LAN_EFFECTIVE_INDEX
        group_index = args.group_index or PAPER_LAN_GROUP_INDEX
    elif args.architecture == "pub-16wdm":
        center_wavelength_nm = args.center_wavelength or DEFAULT_CENTER_WAVELENGTH_NM
        effective_index = args.effective_index or DEFAULT_EFFECTIVE_INDEX
        group_index = args.group_index or DEFAULT_GROUP_INDEX
    else:
        center_wavelength_nm = args.center_wavelength or DEFAULT_CENTER_WAVELENGTH_NM
        effective_index = args.effective_index or DEFAULT_EFFECTIVE_INDEX
        group_index = args.group_index or DEFAULT_GROUP_INDEX

    has_start = args.start_wavelength is not None
    has_stop = args.stop_wavelength is not None
    if has_start != has_stop:
        parser.exit(status=2, message="Both --start-wavelength and --stop-wavelength must be provided together.\n")

    if has_start and has_stop:
        start_nm = args.start_wavelength
        stop_nm = args.stop_wavelength
    else:
        default_span_nm = 110.0 if args.architecture == "pub-16wdm" else 20.0
        span_nm = args.span if args.span is not None else default_span_nm
        start_nm = center_wavelength_nm - span_nm / 2.0
        stop_nm = center_wavelength_nm + span_nm / 2.0

    if args.architecture == "simple":
        design = build_cascaded_mzi(
            center_wavelength_nm=center_wavelength_nm,
            effective_index=effective_index,
            group_index=group_index,
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

    if args.architecture == "pub-16wdm":
        coupler_excess_loss_db = (
            PUB_16WDM_COUPLER_EXCESS_LOSS_DB
            if args.coupler_excess_loss is None
            else args.coupler_excess_loss
        )
        if args.coupler_file is None:
            coupler_response = constant_coupler_response(
                coupling_ratio=args.coupler_ratio,
                excess_loss_db=coupler_excess_loss_db,
                source_name=f"constant k={args.coupler_ratio:.3f}, loss={coupler_excess_loss_db:.3f} dB/coupler",
            )
        else:
            try:
                coupler_response = load_coupler_response_csv(
                    path=args.coupler_file,
                    ratio_column=args.coupler_ratio_column,
                    wavelength_column=args.coupler_wavelength_column,
                    loss_column=args.coupler_loss_column,
                    fallback_loss_db=coupler_excess_loss_db,
                )
            except ValueError as exc:
                parser.exit(status=2, message=f"{exc}\n")

        design = build_pub_16wdm_demux(
            center_wavelength_nm=center_wavelength_nm,
            effective_index=effective_index,
            group_index=group_index,
            channel_spacing_ghz=args.channel_spacing_ghz,
            coupler_response=coupler_response,
            optimize_phase_offsets=not args.no_phase_offset,
        )

        print(format_cascaded_wdm_report(design))
        rows = cascaded_wdm_spectrum(design, start_nm=start_nm, stop_nm=stop_nm, points=args.points)
        peaks = [
            max(rows, key=lambda item, port_index=port_index: item[1][port_index])
            for port_index in range(design.channel_count)
        ]
        print()
        print("Spectrum summary")
        print("-" * 45)
        print(f"Sample range      : {start_nm:.3f} - {stop_nm:.3f} nm")
        for port_index, peak in enumerate(peaks, start=1):
            channel_index = design.port_channels[port_index - 1]
            port_power = peak[1][port_index - 1]
            print(
                f"Port {port_index:>2} / ch {channel_index:>2} : "
                f"{linear_to_db(port_power):>6.2f} dB @ {peak[0]:.3f} nm"
            )

        if args.csv is not None:
            args.csv.parent.mkdir(parents=True, exist_ok=True)
            save_multiport_spectrum_csv(args.csv, rows, port_count=design.channel_count)
            print(f"CSV saved         : {args.csv}")

        if args.plot or args.plot_file is not None:
            try:
                plot_cascaded_wdm_spectrum(
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

    if args.architecture == "lan-wdm":
        design = build_lan_wdm_demux(
            center_wavelength_nm=center_wavelength_nm,
            effective_index=effective_index,
            group_index=group_index,
            channel_spacing_nm=args.channel_spacing,
            coupler_ratio=args.coupler_ratio,
            optimize_phase_offsets=not args.no_phase_offset,
        )

        print(format_lan_demux_report(design))
        rows = lan_demux_spectrum(design, start_nm=start_nm, stop_nm=stop_nm, points=args.points)
        peaks = [
            max(rows, key=lambda item, port_index=port_index: item[1][port_index])
            for port_index in range(8)
        ]
        print()
        print("Spectrum summary")
        print("-" * 45)
        print(f"Sample range      : {start_nm:.3f} - {stop_nm:.3f} nm")
        for port_index, peak in enumerate(peaks, start=1):
            port_power = peak[1][port_index - 1]
            print(f"Port {port_index:>2} peak      : {linear_to_db(port_power):>6.2f} dB @ {peak[0]:.3f} nm")

        if args.csv is not None:
            args.csv.parent.mkdir(parents=True, exist_ok=True)
            save_lan_demux_spectrum_csv(args.csv, rows)
            print(f"CSV saved         : {args.csv}")

        if args.plot or args.plot_file is not None:
            try:
                plot_lan_demux_spectrum(
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
        center_wavelength_nm=center_wavelength_nm,
        effective_index=effective_index,
        group_index=group_index,
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
