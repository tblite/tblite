# This file is part of tblite.
# SPDX-Identifier: LGPL-3.0-or-later
#
# tblite is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tblite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with tblite.  If not, see <https://www.gnu.org/licenses/>.
"""Tests for the qcelemental interface."""

from typing import Any, Dict, Optional, TYPE_CHECKING

import numpy as np
import pytest

try:
    from tblite.qcschema import qcel_v1, qcel_v2, run_schema
except ModuleNotFoundError:
    qcel_v1 = None
    qcel_v2 = None

v1_available = pytest.mark.skipif(
    qcel_v1 is None, reason="QCSchema v1 not available for py314+"
)
v2_available = pytest.mark.skipif(
    qcel_v2 is None, reason="QCSchema v2 not available in current QCElemental"
)

if TYPE_CHECKING:
    Molecule = qcel_v1.Molecule if qcel_v1 is not None else qcel_v2.Molecule
    AtomicInput = qcel_v1.AtomicInput if qcel_v1 is not None else qcel_v2.AtomicInput


@pytest.fixture(
    params=[pytest.param(1, marks=v1_available), pytest.param(2, marks=v2_available)]
)
def qcsk_version(request):
    return request.param


@pytest.fixture
def multiplicity(request) -> int:
    """Multiplicity fixture."""
    return getattr(request, "param", 1)


@pytest.fixture(params=["ala-xab"])
def molecule(request, multiplicity: int, qcsk_version: int) -> "Molecule":
    """Get a molecule for testing."""
    if qcsk_version == 1:
        Molecule = qcel_v1.Molecule
    elif qcsk_version == 2:
        Molecule = qcel_v2.Molecule

    if request.param == "ala-xab":
        return Molecule(
            symbols=list("NHCHCCHHHOCCHHHONHCHHH"),
            geometry=np.array(
                [
                    [+2.65893135608838, -2.39249423371715, -3.66065400053935],
                    [+3.49612941769371, -0.88484673975624, -2.85194146578362],
                    [-0.06354076626069, -2.63180732150005, -3.28819116275323],
                    [-1.07444177498884, -1.92306930149582, -4.93716401361053],
                    [-0.83329925447427, -5.37320588052218, -2.81379718546920],
                    [-0.90691285352090, -1.04371377845950, -1.04918016247507],
                    [-2.86418317801214, -5.46484901686185, -2.49961410229771],
                    [-0.34235262692151, -6.52310417728877, -4.43935278498325],
                    [+0.13208660968384, -6.10946566962768, -1.15032982743173],
                    [-2.96060093623907, +0.01043357425890, -0.99937552379387],
                    [+3.76519127865000, -3.27106236675729, -5.83678272799149],
                    [+6.47957316843231, -2.46911747464509, -6.21176914665408],
                    [+7.32688324906998, -1.67889171278096, -4.51496113512671],
                    [+6.54881843238363, -1.06760660462911, -7.71597456720663],
                    [+7.56369260941896, -4.10015651865148, -6.82588105651977],
                    [+2.64916867837331, -4.60764575400925, -7.35167957128511],
                    [+0.77231592220237, -0.92788783332000, +0.90692539619101],
                    [+2.18437036702702, -2.20200039553542, +0.92105755612696],
                    [+0.01367202674183, +0.22095199845428, +3.27728206652909],
                    [+1.67849497305706, +0.53855308534857, +4.43416031916610],
                    [-0.89254709011762, +2.01704896333243, +2.87780123699499],
                    [-1.32658751691561, -0.95404596601807, +4.30967630773603],
                ]
            ),
            molecular_multiplicity=multiplicity,
        )

    raise ValueError(f"Unknown molecule: {request.param}")


@pytest.fixture(params=["energy", "gradient"])
def driver(request) -> str:
    """Driver fixture."""
    return request.param


@pytest.fixture(params=["GFN1-xTB", "GFN2-xTB"])
def method(request) -> str:
    """Method fixture."""
    return request.param


def get_atomic_input(
    version: int,
    molecule: Dict[str, Any],
    driver: str,
    model: str,
    keywords: Optional[Dict[str, Any]] = None,
    qcel_object: bool = False,
):
    keywords = {} if keywords is None else keywords
    spec = {
        "driver": driver,
        "model": model,
        "keywords": keywords,
    }

    if version == 1:
        input_data = {
            "molecule": molecule,
            **spec,
        }
        if qcel_object:
            return qcel_v1.AtomicInput(**input_data)
        return input_data

    if version == 2:
        input_data = {
            "molecule": molecule,
            "specification": spec,
        }
        if qcel_object:
            return qcel_v2.AtomicInput(**input_data)
        return input_data

    raise ValueError(f"Unsupported version: {version}")


@pytest.fixture()
def atomic_input(
    molecule: "Molecule", driver: str, method: str, qcsk_version: int
) -> "AtomicInput":
    """AtomicInput fixture."""
    return get_atomic_input(
        qcsk_version,
        molecule=molecule,
        driver=driver,
        model={"method": method},
        keywords={"spin-polarization": 1.0},
        qcel_object=True,
    )


@pytest.fixture()
def return_result(molecule: "Molecule", driver: str, method: str) -> Any:
    """Return result fixture."""
    if qcel_v1 is None and qcel_v2 is None:
        return None

    # fmt: off
    return {
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "energy",
            "GFN1-xTB",
        ): -34.98079481580359,
        (
            "65a3bab7309579268e383249bad8d33095c32db5", 
            "energy",
            "GFN1-xTB",
        ): -34.79969349420318,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "energy",
            "GFN2-xTB",
        ): -32.96247200752864,
        (
            "65a3bab7309579268e383249bad8d33095c32db5", 
            "energy",
            "GFN2-xTB",
        ): -32.79550096196652,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "gradient",
            "GFN1-xTB",
        ): np.array(
            [
                [-7.8113833723432776e-3, 1.1927359496178194e-3, 4.5384293534289468e-3],
                [4.6431896826996466e-4, -1.1893457514353986e-3, -2.0196769165192010e-3],
                [1.4651707521730923e-3, -2.1095273552157868e-3, -1.9374661833396272e-3],
                [1.1712050542529663e-3, -6.2696965846751425e-4, 6.1050802435153369e-3],
                [-1.1270171972379136e-3, 9.4432659541408457e-4, -2.0020105757037662e-3],
                [1.1646562626373442e-2, -7.8838363893641971e-3, -9.4098748734209019e-3],
                [4.1886781466182168e-4, -2.0573790837090546e-4, 4.0632563116798048e-4],
                [-1.8053991102095574e-4, 1.1681264331977676e-3, 1.6828801638121846e-3],
                [-5.9900194915516749e-4, 3.4776846711688584e-4, -3.7491749453125091e-4],
                [-1.4520319022330924e-2, 6.0131467543009937e-3, 7.5687331646021375e-4],
                [1.5602405715421146e-2, 7.9473817513023640e-3, 8.0623888147124366e-3],
                [-3.2563959623975588e-4, -2.1680106787522012e-4, -8.8683653969162549e-4],
                [7.5180811527270801e-4, -1.9128778304917517e-4, -1.0174498970762392e-3],
                [7.4132920268697234e-4, -6.3819759962106609e-4, 5.1853393972177886e-4],
                [-7.5646645647444864e-4, 1.5490223231011606e-3, 6.5407919650053525e-4],
                [-9.0634683701016835e-3, -8.8383472527482337e-3, -1.2123112366846918e-2],
                [3.1541524835559096e-3, 2.7233221491356533e-3, 1.0127030544629243e-2],
                [-1.4266036482687263e-3, 1.2132816331079002e-3, -2.7113843360055362e-3],
                [2.0255265706568555e-4, -1.3134487584798992e-3, -7.9291928858605555e-4],
                [-2.3317056624669709e-3, -2.1233032658492240e-4, -1.1563077745307797e-3],
                [1.6004440124424719e-3, -3.2444847372739881e-3, 1.6000422202203360e-3],
                [9.2332778346354430e-4, 3.5712025321916925e-3, -1.9707177917017700e-5],
            ]
        ),
        (
            "65a3bab7309579268e383249bad8d33095c32db5",
            "gradient",
            "GFN1-xTB",
        ): np.array(
            [
                [5.8524614928954187e-3, -2.8273860862564008e-3, -2.5050614815354293e-2],
                [1.1339401892754533e-3, -3.8043823659953784e-3, 4.6415780745915643e-4],
                [6.1161622297322229e-4, 4.2685677262145104e-3, 1.4570851927090144e-3],
                [3.8586597902057181e-3, -2.5177483893864799e-3, 4.5515645848229879e-3],
                [-2.8918206405656940e-3, 2.8455866696309368e-3, -5.9668370428631338e-4],
                [-3.1705370345825074e-2, 2.4083262454045112e-2, 2.7384582155083294e-2],
                [-5.1707111952240632e-4, -2.5643418411340111e-4, 1.2581018565013965e-3],
                [-1.2470352203053496e-3, 2.5614875177244435e-3, 4.3056402494263653e-3],
                [-6.6844836248905481e-4, 2.8089761221454008e-4, 7.8809940986204145e-4],
                [6.2706334682231246e-2, -3.4003886955687332e-2, 8.0586129937534315e-4],
                [-9.4087011001570357e-3, -3.4501995852940378e-3, -1.7483270978818424e-4],
                [-9.2207770819250441e-4, -3.5421985555735337e-3, 9.5662473713666832e-4],
                [1.3078320300956743e-3, 4.1263243281132945e-5, -1.3500611258690431e-3],
                [-9.0178431942333195e-4, -2.0509954569158484e-3, 1.5334520471585486e-3],
                [-6.9184533697939866e-4, 2.1987490376191906e-3, 4.6577555215977903e-4],
                [9.1755875101464449e-3, 1.5842405697885673e-2, 1.4515316734675897e-2],
                [-3.8904615489211355e-2, -1.1230078292625741e-2, -3.0192692178494903e-2],
                [3.2109705087129807e-4, 5.0278073973231176e-3, -4.7478466044375703e-3],
                [1.6316664176902350e-3, 5.0025666147582287e-3, 3.8013881083669714e-3],
                [-2.8338123665745895e-3, 1.1958861407827702e-3, -6.8946026182785858e-4],
                [2.3970530203739990e-3, -3.0498463357009459e-3, 2.4297750760109019e-3],
                [1.6963336024870001e-3, 3.3846760960695773e-3, -1.9152334106901079e-3],
            ]
        ),
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "gradient",
            "GFN2-xTB",
        ): np.array(
            [
                [-2.8793799913136380e-3, 3.6460501266212778e-3, 1.0760791938638532e-2],
                [4.4678599343911276e-4, -4.1389479945372787e-4, -1.6588444399860403e-3],
                [-5.2514562332919819e-3, -3.5790236939509104e-3, -3.7701781817474711e-3],
                [1.6246105310470146e-3, 3.3329849316788402e-4, 5.2523889343070009e-3],
                [2.7987398962996536e-4, 1.4306153660053899e-3, -2.4860802541022639e-3],
                [1.3079422937187748e-2, -6.5673919526485743e-3, -9.0220988900455105e-3],
                [-1.1471778995347530e-3, 5.7911970791709148e-5, 5.9833508963369865e-4],
                [-1.2272095567273031e-4, -3.1971719956945088e-4, 5.0170686443738396e-4],
                [1.4364574920106048e-4, 3.3059193829213012e-4, 7.6579452785142488e-4],
                [-1.5793744333963813e-2, 5.5330915185388980e-3, 1.7671683785165003e-4],
                [1.2223295427589707e-2, 6.4327672083176893e-3, 7.2755941232891530e-3],
                [-9.0551792272509739e-4, 2.4507458199862330e-5, -1.5354579310736880e-3],
                [1.0467950926876190e-3, 3.9864080026996094e-4, -1.2505277667061539e-4],
                [2.8780762038940164e-4, 8.3414968003247950e-5, 2.8920334379877398e-5],
                [-6.5409450906848499e-4, 1.9393428922497808e-4, 6.3099484385291252e-4],
                [-8.1386141087597102e-3, -8.7642669443246026e-3, -1.2729326807486986e-2],
                [1.0141941984778097e-2, -1.1848858357028618e-3, 5.2114368007609652e-3],
                [-9.1262500476170961e-4, 5.8341311743242258e-4, -2.0431514469271190e-3],
                [-4.4227155696016868e-3, 1.5005306210378255e-3, 1.5457167117877691e-3],
                [-1.1435610900201417e-3, -2.1401635839282676e-4, -6.6270887861705056e-4],
                [1.3141203673661993e-3, -1.3944503339669997e-3, 1.6665679713092187e-3],
                [7.8330792539781752e-4, 1.8888792421067396e-3, -3.8206537144275390e-4],
            ],
        ),
        (
            "65a3bab7309579268e383249bad8d33095c32db5", 
            "gradient",
            "GFN2-xTB",
        ): np.array(
            [
                [8.9047472041766714e-3, 9.5607884444064183e-4, -2.0207382561384348e-2],
                [1.2421997624991524e-3, -4.1923149859461840e-3, 5.7191874243116270e-4],
                [-2.2574968847484576e-3, 1.6557395021578063e-3, 8.0267045689306820e-4],
                [5.2249531945359811e-3, -1.1939260781055974e-3, 4.8042424956936929e-3],
                [-1.9782679250485149e-3, 1.3446000662190099e-3, -2.2001974385799923e-3],
                [-3.9738793437907795e-2, 2.1197953201955746e-2, 1.5247311907206403e-2],
                [-1.2242960187947738e-3, 2.2697432798506555e-4, 1.5636237394129760e-3],
                [-1.5355554337369052e-3, 1.1773534096960524e-3, 3.9239318822003773e-3],
                [-2.5685691169087657e-4, 6.5672081371170182e-4, 1.4060135522050222e-3],
                [5.7184838051717493e-2, -2.9817033711430491e-2, 3.5460781848660520e-3],
                [-1.6884748580541110e-2, -1.2555342101938205e-2, -1.2838662695503785e-2],
                [-1.2311549320041084e-3, -3.0607127282149742e-3, 1.5888644497996985e-3],
                [2.3019232268532133e-3, 9.7165262841857351e-4, -4.5119341880787332e-4],
                [-1.5506098550687828e-3, -1.9369010122701987e-3, 1.7871259175595224e-3],
                [-5.8808020643440532e-4, 1.1369206741758716e-3, 3.9542398353661321e-4],
                [1.4834642395475424e-2, 2.3171870918960714e-2, 2.4016664981649165e-2],
                [-2.3690908686237466e-2, -1.1167134989771683e-2, -2.4170736382859302e-2],
                [1.7138999209587496e-3, 3.5315927824726831e-3, -3.7555205492248679e-3],
                [-2.4669122256482951e-3, 6.9132394242477774e-3, 4.8203547434331109e-3],
                [-1.4685833418781526e-3, 7.8984853258663224e-4, -1.8315220150681287e-4],
                [1.5859948906085583e-3, -1.6802967562846727e-3, 2.2718006537276250e-3],
                [1.8790657929145008e-3, 1.8731172369336867e-3, -2.9391804427475331e-3],
            ],
        ),
    }[(molecule.get_hash(), driver, method)]
    # fmt: on


@pytest.mark.skipif(qcel_v1 is None and qcel_v2 is None, reason="requires qcelemental")
@pytest.mark.parametrize("multiplicity", [1, 3], indirect=True)
def test_qcschema(atomic_input: "AtomicInput", return_result: Any) -> None:
    """Test qcschema interface."""
    atomic_result = run_schema(atomic_input)

    assert atomic_result.success
    assert pytest.approx(atomic_result.return_result) == return_result


@pytest.fixture(
    params=[
        {"cosmo-solvation": "methanol"},
        {"cosmo-solvation": 7.0},
        {"cpcm-solvation": 7.0},
        {"pcm-solvation": 7.0},
        {"alpb-solvation": ["water", "bar1mol"]},
        {"gbsa-solvation": ["methanol", "reference"]},
        {"gbe-solvation": [7.0, "p16"]},
        {"gb-solvation": [7.0, "still"]},
    ]
)
def solvation(request) -> dict:
    """Solvation fixture."""
    return request.param


@pytest.fixture()
def atomic_input_solvation(
    molecule: "Molecule", method: str, solvation: dict, qcsk_version: int
) -> "AtomicInput":
    """AtomicInput fixture."""
    return get_atomic_input(
        qcsk_version,
        molecule=molecule,
        driver="energy",
        model={"method": method},
        keywords=solvation,
        qcel_object=False,
    )


@pytest.fixture()
def return_result_solvation(molecule: "Molecule", method: str, solvation: dict) -> Any:
    """Return result fixture."""
    if qcel_v1 is None and qcel_v2 is None:
        return None

    # fmt: off
    return {
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "cosmo-solvation",
            "methanol",
            "GFN1-xTB",
        ): -35.00448615440448,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "cosmo-solvation",
            "methanol",
            "GFN2-xTB",
        ): -32.98232152157144,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "cosmo-solvation",
            7.0,
            "GFN1-xTB",
        ): -34.999904668295585,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "cosmo-solvation",
            7.0,
            "GFN2-xTB",
        ): -32.97812262050184,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "cpcm-solvation",
            7.0,
            "GFN1-xTB",
        ): -35.00155428273607,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "cpcm-solvation",
            7.0,
            "GFN2-xTB",
        ): -32.979610392987865,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "pcm-solvation",
            7.0,
            "GFN1-xTB",
        ): -34.99927517635751,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "pcm-solvation",
            7.0,
            "GFN2-xTB",
        ): -32.97745476937624,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "alpb-solvation",
            ("water", "bar1mol"),
            "GFN1-xTB",
        ): -34.9966872968042,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "alpb-solvation",
            ("water", "bar1mol"),
            "GFN2-xTB",
        ): -32.9777888674693,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "gbsa-solvation",
            ("methanol", "reference"),
            "GFN1-xTB",
        ): -34.9989058302890,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "gbsa-solvation",
            ("methanol", "reference"),
            "GFN2-xTB",
        ): -32.9763255038126,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "gbe-solvation",
            (7.0, "p16"),
            "GFN1-xTB",
        ): -34.9859026499920,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "gbe-solvation",
            (7.0, "p16"),
            "GFN2-xTB",
        ): -32.9667327170214,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "gb-solvation",
            (7.0, "still"),
            "GFN1-xTB",
        ): -34.9859968024425,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "gb-solvation",
            (7.0, "still"),
            "GFN2-xTB",
        ): -32.9668850575386,
    }[
        (
            molecule.get_hash(),
            list(solvation.keys())[0],
            tuple(list(solvation.values())[0])
            if isinstance(list(solvation.values())[0], list)
            else list(solvation.values())[0],
            method,
        )
    ]
    # fmt: on


@pytest.mark.skipif(qcel_v1 is None and qcel_v2 is None, reason="requires qcelemental")
def test_qcschema_solvation(
    atomic_input_solvation: "AtomicInput", return_result_solvation: Any
) -> None:
    """Test qcschema interface."""
    atomic_result = run_schema(atomic_input_solvation)

    assert atomic_result.success
    assert pytest.approx(atomic_result.return_result) == return_result_solvation


@pytest.mark.skipif(qcel_v1 is None and qcel_v2 is None, reason="requires qcelemental")
def test_unsupported_driver(molecule: "Molecule", qcsk_version: int):
    """Test unsupported driver name."""
    atomic_inp = get_atomic_input(
        qcsk_version,
        molecule=molecule,
        driver="hessian",
        model={"method": "GFN1-xTB"},
    )

    atomic_result = run_schema(atomic_inp)

    assert not atomic_result.success
    assert atomic_result.error.error_type == "input_error"
    assert (
        "Driver 'hessian' is not supported by tblite."
        in atomic_result.error.error_message
    )


@pytest.mark.skipif(qcel_v1 is None and qcel_v2 is None, reason="requires qcelemental")
def test_unsupported_method(molecule: "Molecule", qcsk_version: int):
    """Test unsupported method name."""
    atomic_inp = get_atomic_input(
        qcsk_version,
        molecule=molecule,
        driver="energy",
        model={"method": "GFN-xTB"},
    )

    atomic_result = run_schema(atomic_inp)

    assert not atomic_result.success
    assert atomic_result.error.error_type == "input_error"
    assert (
        "Model 'GFN-xTB' is not supported by tblite."
        in atomic_result.error.error_message
    )


@pytest.mark.skipif(qcel_v1 is None and qcel_v2 is None, reason="requires qcelemental")
def test_unsupported_basis(molecule: "Molecule", qcsk_version: int):
    """Test unsupported basis set."""
    atomic_inp = get_atomic_input(
        qcsk_version,
        molecule=molecule,
        driver="energy",
        model={"method": "GFN1-xTB", "basis": "def2-SVP"},
    )

    atomic_result = run_schema(atomic_inp)

    assert not atomic_result.success
    assert atomic_result.error.error_type == "input_error"
    assert (
        "Basis sets are not supported by tblite." in atomic_result.error.error_message
    )


@pytest.mark.skipif(qcel_v1 is None and qcel_v2 is None, reason="requires qcelemental")
def test_unsupported_keywords(molecule: "Molecule", qcsk_version: int):
    """Test unsupported keywords."""
    atomic_inp = get_atomic_input(
        qcsk_version,
        molecule=molecule,
        driver="gradient",
        model={"method": "GFN1-xTB"},
        keywords={"unsupported": True},
    )

    atomic_result = run_schema(atomic_inp)

    assert not atomic_result.success
    assert atomic_result.error.error_type == "input_error"
    assert "Unknown keywords: unsupported" in atomic_result.error.error_message


@pytest.mark.skipif(qcel_v1 is None and qcel_v2 is None, reason="requires qcelemental")
def test_scf_not_converged(molecule: "Molecule", qcsk_version: int):
    """Test unconverged SCF."""
    atomic_inp = get_atomic_input(
        qcsk_version,
        molecule=molecule,
        driver="gradient",
        model={"method": "GFN1-xTB"},
        keywords={"max-iter": 3},
    )

    atomic_result = run_schema(atomic_inp)

    assert not atomic_result.success
    assert atomic_result.error.error_type == "execution_error"
    assert "SCF not converged in 3 cycles" in atomic_result.error.error_message
