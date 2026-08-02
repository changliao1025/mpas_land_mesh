import importlib
import sys


def test_gdal_imports_gracefully_when_native_extension_is_unavailable(monkeypatch):
    monkeypatch.setitem(sys.modules, 'osgeo', None)
    monkeypatch.setitem(sys.modules, 'osgeo.gdal', None)
    monkeypatch.setitem(sys.modules, 'osgeo.ogr', None)
    monkeypatch.setitem(sys.modules, 'osgeo.osr', None)

    import mpas_land_mesh.utilities.vector as vector
    assert vector.gdal is None
    assert vector.ogr is None
    assert vector.osr is None
