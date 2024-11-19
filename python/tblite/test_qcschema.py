from typing import Any

import numpy as np
import pytest

try:
    import qcelemental as qcel
    from qcelemental.models import AtomicInput, Molecule
    from tblite.qcschema import run_schema
except ModuleNotFoundError:
    qcel = None
    AtomicInput = dict
    Molecule = dict


@pytest.fixture(params=["ala-xab"])
def molecule(request) -> Molecule:
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
        )

    raise ValueError(f"Unknown molecule: {request.param}")


@pytest.fixture(params=["energy", "gradient"])
def driver(request) -> str:
    return request.param


@pytest.fixture(params=["GFN1-xTB", "GFN2-xTB"])
def method(request) -> str:
    return request.param


@pytest.fixture()
def atomic_input(molecule: Molecule, driver: str, method: str) -> AtomicInput:
    return AtomicInput(
        molecule=molecule,
        driver=driver,
        model={"method": method},
    )


@pytest.fixture()
def return_result(molecule: Molecule, driver: str, method: str) -> Any:
    if qcel is None:
        return None
    return {
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "energy",
            "GFN1-xTB",
        ): -34.98079481580359,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "energy",
            "GFN2-xTB",
        ): -32.96247200752864,
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
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "gradient",
            "GFN2-xTB",
        ): np.array(
            [
                [-2.8763728770498059e-3, 3.6434117451340639e-3, 1.0752016108561982e-2],
                [4.4678199702658915e-4, -4.1415548977460947e-4, -1.6583158023285350e-3],
                [-5.2507431234969960e-3, -3.5785920414832134e-3, -3.7696789673926740e-3],
                [1.6244739083799741e-3, 3.3338843276814147e-4, 5.2523752797833641e-3],
                [2.7999803310429892e-4, 1.4312674725036575e-3, -2.4858112663419798e-3],
                [1.3090058013057956e-2, -6.5715625689103368e-3, -9.0191025142361859e-3],
                [-1.1473355542259416e-3, 5.7814760290257515e-5, 5.9832857670227466e-4],
                [-1.2273010075096477e-4, -3.1980940288597041e-4, 5.0156274813269695e-4],
                [1.4368065055278624e-4, 3.3050603548747309e-4, 7.6584802762956364e-4],
                [-1.5800077393544824e-2, 5.5365186088464999e-3, 1.7682778067956641e-4],
                [1.2225203364730703e-2, 6.4403779555781172e-3, 7.2872510746961267e-3],
                [-9.0925674185912101e-4, 2.3218888809526880e-5, -1.5348977326777743e-3],
                [1.0468938953083059e-3, 3.9870261221090688e-4, -1.2509198028748229e-4],
                [2.8775846835979148e-4, 8.3392790833921953e-5, 2.8790511400042136e-5],
                [-6.5400300154964302e-4, 1.9383201778978590e-4, 6.3086666998487021e-4],
                [-8.1410643777397960e-3, -8.7674163235359689e-3, -1.2732857726175186e-2],
                [1.0137767740747086e-2, -1.1853820088445544e-3, 5.2080908879111347e-3],
                [-9.1245045379595587e-4, 5.8368999845378263e-4, -2.0429521928185526e-3],
                [-4.4224663957623629e-3, 1.5003741790590788e-3, 1.5449104035364322e-3],
                [-1.1434980317896989e-3, -2.1396289672603825e-4, -6.6261217226958405e-4],
                [1.3140684680651973e-3, -1.3944469817754714e-3, 1.6666081386766634e-3],
                [7.8331351223255292e-4, 1.8888322161710286e-3, -3.8215585316686302e-4],
            ],
        ),
    }[(molecule.get_hash(), driver, method)]


@pytest.mark.skipif(qcel is None, reason="requires qcelemental")
def test_qcschema(atomic_input: AtomicInput, return_result: Any) -> None:
    atomic_result = run_schema(atomic_input)

    assert atomic_result.success
    print(atomic_result.return_result)
    assert pytest.approx(atomic_result.return_result) == return_result


@pytest.mark.skipif(qcel is None, reason="requires qcelemental")
def test_unsupported_driver(molecule: Molecule):
    atomic_input = AtomicInput(
        molecule=molecule,
        driver="hessian",
        model={"method": "GFN1-xTB"},
    )

    atomic_result = run_schema(atomic_input)

    assert not atomic_result.success
    assert atomic_result.error.error_type == "input_error"
    assert "Driver 'hessian' is not supported by tblite." in atomic_result.error.error_message


@pytest.mark.skipif(qcel is None, reason="requires qcelemental")
def test_unsupported_method(molecule: Molecule):
    atomic_input = AtomicInput(
        molecule=molecule,
        driver="energy",
        model={"method": "GFN-xTB"},
    )

    atomic_result = run_schema(atomic_input)

    assert not atomic_result.success
    assert atomic_result.error.error_type == "input_error"
    assert "Model 'GFN-xTB' is not supported by tblite." in atomic_result.error.error_message


@pytest.mark.skipif(qcel is None, reason="requires qcelemental")
def test_unsupported_basis(molecule: Molecule):
    atomic_input = AtomicInput(
        molecule=molecule,
        driver="energy",
        model={"method": "GFN1-xTB", "basis": "def2-SVP"},
    )

    atomic_result = run_schema(atomic_input)

    assert not atomic_result.success
    assert atomic_result.error.error_type == "input_error"
    assert "Basis sets are not supported by tblite." in atomic_result.error.error_message


@pytest.mark.skipif(qcel is None, reason="requires qcelemental")
def test_unsupported_keywords(molecule: Molecule):
    atomic_input = AtomicInput(
        molecule=molecule,
        driver="gradient",
        model={"method": "GFN1-xTB"},
        keywords={"unsupported": True},
    )

    atomic_result = run_schema(atomic_input)

    assert not atomic_result.success
    assert atomic_result.error.error_type == "input_error"
    assert "Unknown keywords: unsupported" in atomic_result.error.error_message


@pytest.mark.skipif(qcel is None, reason="requires qcelemental")
def test_scf_not_converged(molecule: Molecule):
    atomic_input = AtomicInput(
        molecule=molecule,
        driver="gradient",
        model={"method": "GFN1-xTB"},
        keywords={"max-iter": 3},
    )

    atomic_result = run_schema(atomic_input)

    assert not atomic_result.success
    assert atomic_result.error.error_type == "execution_error"
    assert "SCF not converged in 3 cycles" in atomic_result.error.error_message
