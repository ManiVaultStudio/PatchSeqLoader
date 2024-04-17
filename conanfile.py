from conans import ConanFile
from conan.tools.cmake import CMakeDeps, CMake, CMakeToolchain
from conans.tools import save, load
import os
import shutil
import pathlib
import subprocess
from rules_support import PluginBranchInfo


class PatchSeqDataLoaderConan(ConanFile):
    """Class to package the plugin using conan

    Packages both RELEASE and DEBUG.
    Uses rules_support (github.com/hdps/rulessupport) to derive
    versioninfo based on the branch naming convention
    as described in https://github.com/ManiVaultStudio/core/wiki/Branch-naming-rules
    """

    name = "PatchSeqDataLoader"
    description = """Loader for multi-modal patch-seq data, consisting of a gene expression matrix, electrophysiology and morphology data."""
    topics = ("manivault", "plugin", "loader", "patch-seq")
    url = "https://github.com/ManiVaultStudio/PatchSeqLoader"
    author = "julianthijssen@gmail.com"  # conan recipe author
    license = "LGPL 3.0"

    short_paths = True
    generators = "CMakeDeps"

    # Options may need to change depending on the packaged library
    settings = {"os": None, "build_type": None, "compiler": None, "arch": None}
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {"shared": True, "fPIC": True}

    requires = ("CellMorphologyData/1.0@lkeb/stable")

    # Qt requirement is inherited from mv-core

    scm = {
        "type": "git",
        "subfolder": "hdps/PatchSeqDataLoader",
        "url": "auto",
        "revision": "auto",
    }

    def __get_git_path(self):
        path = load(
            pathlib.Path(pathlib.Path(__file__).parent.resolve(), "__gitpath.txt")
        )
        print(f"git info from {path}")
        return path

    def export(self):
        print("In export")
        # save the original source path to the directory used to build the package
        save(
            pathlib.Path(self.export_folder, "__gitpath.txt"),
            str(pathlib.Path(__file__).parent.resolve()),
        )

    def set_version(self):
        # Assign a version from the branch name
        branch_info = PluginBranchInfo(self.recipe_folder)
        self.version = branch_info.version
        # print(f"Got version: {self.version}")

    def requirements(self):
        branch_info = PluginBranchInfo(self.__get_git_path())
        print(f"Core requirement {branch_info.core_requirement}")
        self.requires(branch_info.core_requirement)

    def configure(self):
        pass

    def system_requirements(self):
        #  May be needed for macOS or Linux
        pass

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def generate(self):
        generator = None
        if self.settings.os == "Macos":
            generator = "Xcode"
        if self.settings.os == "Linux":
            generator = "Ninja Multi-Config"
        # Use the Qt provided .cmake files
        qtpath = pathlib.Path(self.deps_cpp_info["qt"].rootpath)
        qt_root = str(list(qtpath.glob("**/Qt6Config.cmake"))[0].parents[3].as_posix())

        tc = CMakeToolchain(self, generator=generator)
        if self.settings.os == "Windows" and self.options.shared:
            tc.variables["CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS"] = True
        if self.settings.os == "Linux" or self.settings.os == "Macos":
            tc.variables["CMAKE_CXX_STANDARD_REQUIRED"] = "ON"
        tc.variables["CMAKE_PREFIX_PATH"] = qt_root
        
        # Set the installation directory for ManiVault based on the MV_INSTALL_DIR environment variable
        # or if none is specified, set it to the build/install dir.
        if not os.environ.get("MV_INSTALL_DIR", None):
            os.environ["MV_INSTALL_DIR"] = os.path.join(self.build_folder, "install")
        print("MV_INSTALL_DIR: ", os.environ["MV_INSTALL_DIR"])
        self.install_dir = pathlib.Path(os.environ["MV_INSTALL_DIR"]).as_posix()
        # Give the installation directory to CMake
        MV_CMD_PATH = pathlib.Path(self.deps_cpp_info["CellMorphologyData"].rootpath).as_posix()
        print(f"MV_CMD_INSTALL_DIR: {MV_CMD_PATH}")
        tc.variables["MV_INSTALL_DIR"] = self.install_dir
        tc.variables["MV_CMD_INSTALL_DIR"] = MV_CMD_PATH
        
        tc.generate()

    def _configure_cmake(self):
        cmake = CMake(self)
        cmake.configure(build_script_folder="hdps/PatchSeqDataLoader")
        cmake.verbose = True
        return cmake

    def build(self):
        print("Build OS is : ", self.settings.os)

        # The BinNIO plugins expect the HDPS package to be in this install dir
        hdps_pkg_root = self.deps_cpp_info["hdps-core"].rootpath
        print("Install dir type: ", self.install_dir)
        shutil.copytree(hdps_pkg_root, self.install_dir)

        cmake = self._configure_cmake()
        cmake.build(build_type="Debug")
        cmake.install(build_type="Debug")

        # cmake_release = self._configure_cmake()
        cmake.build(build_type="Release")
        cmake.install(build_type="Release")

    def package(self):
        package_dir = os.path.join(self.build_folder, "package")
        print("Packaging install dir: ", package_dir)
        subprocess.run(
            [
                "cmake",
                "--install",
                self.build_folder,
                "--config",
                "Debug",
                "--prefix",
                os.path.join(package_dir, "Debug"),
            ]
        )
        subprocess.run(
            [
                "cmake",
                "--install",
                self.build_folder,
                "--config",
                "Release",
                "--prefix",
                os.path.join(package_dir, "Release"),
            ]
        )
        self.copy(pattern="*", src=package_dir)
        # Add the debug support files to the package
        # (*.pdb) if building the Visual Studio version
        if self.settings.compiler == "Visual Studio":
            self.copy("*.pdb", dst="Debug/Plugins", keep_path=False)

    def package_info(self):
        self.cpp_info.debug.libdirs = ["Debug/lib"]
        self.cpp_info.debug.bindirs = ["Debug/Plugins", "Debug"]
        self.cpp_info.debug.includedirs = ["Debug/include", "Debug"]
        self.cpp_info.release.libdirs = ["Release/lib"]
        self.cpp_info.release.bindirs = ["Release/Plugins", "Release"]
        self.cpp_info.release.includedirs = ["Release/include", "Release"]
