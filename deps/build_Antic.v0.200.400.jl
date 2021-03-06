using BinaryProvider # requires BinaryProvider 0.3.0 or later

# Parse some basic command-line arguments
const verbose = "--verbose" in ARGS
const prefix = Prefix(get([a for a in ARGS if a != "--verbose"], 1, joinpath(@__DIR__, "usr")))
products = [
    LibraryProduct(prefix, ["libantic"], :libantic),
]

# Download binaries from hosted location
bin_prefix = "https://github.com/JuliaBinaryWrappers/Antic_jll.jl/releases/download/Antic-v0.200.400+0"

# Listing of files generated by BinaryBuilder:
download_info = Dict(
    Linux(:aarch64, libc=:glibc) => ("$bin_prefix/Antic.v0.200.400.aarch64-linux-gnu.tar.gz", "6ddd116dc0dbc19094f84ba3244094d927990b1872d04a352fc4a1eb2867a97d"),
    Linux(:aarch64, libc=:musl) => ("$bin_prefix/Antic.v0.200.400.aarch64-linux-musl.tar.gz", "6f4878572d8695c99980049c3f0c9037a44b2e1687318f493166cdbcff664408"),
    Linux(:armv7l, libc=:glibc, call_abi=:eabihf) => ("$bin_prefix/Antic.v0.200.400.armv7l-linux-gnueabihf.tar.gz", "e5cfc54cc8f9549c821f597c41318a7fa92c3e1daf3c8fe6f27e85eb0172249a"),
    Linux(:armv7l, libc=:musl, call_abi=:eabihf) => ("$bin_prefix/Antic.v0.200.400.armv7l-linux-musleabihf.tar.gz", "9408fa4516ced962a86a38e65777e6c18a2be7ec053a7ee9095c8d1c9c87ba8d"),
    Linux(:i686, libc=:glibc) => ("$bin_prefix/Antic.v0.200.400.i686-linux-gnu.tar.gz", "2987d667b01d8ab969b6007937f53ff52fdc5e3d0a93a96e844d8964f0300dab"),
    Linux(:i686, libc=:musl) => ("$bin_prefix/Antic.v0.200.400.i686-linux-musl.tar.gz", "8932d17b82315ead981214ce8142542fdf07283b497efd0cefaf398984269091"),
    Windows(:i686) => ("$bin_prefix/Antic.v0.200.400.i686-w64-mingw32.tar.gz", "4df374e094e7837a4e520357d30512c54a3318bcadfe21a102677a625989bfb0"),
    Linux(:powerpc64le, libc=:glibc) => ("$bin_prefix/Antic.v0.200.400.powerpc64le-linux-gnu.tar.gz", "fb05c69f0f34c2b115e38dbe96f526a8cde8b57cbde46309e9c56a813a23a347"),
    MacOS(:x86_64) => ("$bin_prefix/Antic.v0.200.400.x86_64-apple-darwin.tar.gz", "ddc2ff2b2e4cfe5e160bfa6f4edcb1025227d102946feea17eeb3867cb03af59"),
    Linux(:x86_64, libc=:glibc) => ("$bin_prefix/Antic.v0.200.400.x86_64-linux-gnu.tar.gz", "93fae5f365f019630fd24dc021b9257d301d80773810edde788d63456ca762c1"),
    Linux(:x86_64, libc=:musl) => ("$bin_prefix/Antic.v0.200.400.x86_64-linux-musl.tar.gz", "33d77051af64ec1bde960acf5935889193955060f496d3753eccdc892f36efb3"),
    FreeBSD(:x86_64) => ("$bin_prefix/Antic.v0.200.400.x86_64-unknown-freebsd.tar.gz", "d4cf8db555fe466550e0e924c9c77f8437301dfc1d6f7c3083e2f0e299461c74"),
    Windows(:x86_64) => ("$bin_prefix/Antic.v0.200.400.x86_64-w64-mingw32.tar.gz", "3a3b5aab74c300446bae6598a6ae0719d78f6dd028dc1f13a1888d687a973029"),
)

# Install unsatisfied or updated dependencies:
unsatisfied = any(!satisfied(p; verbose=verbose) for p in products)
dl_info = choose_download(download_info, platform_key_abi())
if dl_info === nothing && unsatisfied
    # If we don't have a compatible .tar.gz to download, complain.
    # Alternatively, you could attempt to install from a separate provider,
    # build from source or something even more ambitious here.
    error("Your platform (\"$(Sys.MACHINE)\", parsed as \"$(triplet(platform_key_abi()))\") is not supported by this package!")
end

# If we have a download, and we are unsatisfied (or the version we're
# trying to install is not itself installed) then load it up!
if unsatisfied || !isinstalled(dl_info...; prefix=prefix)
    # Download and install binaries
    install(dl_info...; prefix=prefix, force=true, verbose=verbose)
end

# Write out a deps.jl file that will contain mappings for our products
write_deps_file(joinpath(@__DIR__, "deps.jl"), products, verbose=verbose)
