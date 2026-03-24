plugins {
    id("com.gradleup.shadow") version "8.3.5"
    id("qupath-conventions")
}

qupathExtension {
    name = "FlowPath"
    group = "io.github.qupath"
    version = "1.0.0"
    description = "Interactive tree-based cell phenotyping gating for multiplexed imaging"
    automaticModule = "qupath.ext.flowpath"
}

dependencies {
    shadow(libs.bundles.qupath)
    shadow(libs.bundles.logging)
    shadow(libs.qupath.fxtras)
    testImplementation(libs.bundles.qupath)
    testImplementation(libs.junit)
    testRuntimeOnly("org.junit.platform:junit-platform-launcher")
}

tasks.test {
    useJUnitPlatform()
}
